//!
//! Copied and modified from rmpfit crate Copyright (c) Vadim Dyadkin
//! Rust implementation of [CMPFIT](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html)
//!
use crate::fit::config::MPFitConfig;
use crate::fit::jacobian::JacMatrix;
use crate::fit::status::MPFitStatus;
use crate::fit::ParameterBounds;
use num::traits::NumAssign;

use crate::constants::FloatConst;
use crate::fit::enorm::ENorm;
use crate::fit::enums::{MPFitError, MPFitInfo};
use crate::fit::ImpedanceModel;

pub struct MPFit<'a, T, U> {
    /// Number of residuals (2 * number of frequency points: real and
    /// imaginary parts of the weighted impedance deviation are stacked).
    m: usize,
    /// Number of parameters.
    npar: usize,
    /// Number of free parameters.
    nfree: usize,
    /// Indices into `xall` of the free (non-fixed) parameters; length `nfree`.
    ifree: Vec<usize>,
    /// Array of length `m` which contains the functions evaluated at the output `x`.
    fvec: Vec<T>,
    /// Variable set to the number of calls to the model function.
    nfev: usize,
    /// Full parameter vector, length `npar`. Scratch buffer used while
    /// the Jacobian and trial steps are being computed.
    xnew: Vec<T>,
    /// Vector of free parameter values, length `nfree`. This is what
    /// the LM step actually operates on — the trust region, Jacobian columns,
    /// and damping all live in this `nfree`-dimensional space. Maps back to
    /// the user-facing layout via `xall[ifree[j]] = x[j]`.
    x: Vec<T>,
    /// Full parameter vector, length `npar`. Holds both fixed and free
    /// parameters and is what gets handed to `model.evaluate`.
    /// During a step, the in-progress full vector lives in
    /// `xnew` while `xall` only reflects accepted values.
    xall: Vec<T>,
    /// Array of length n which contains the first n elements of the
    /// vector `(q transpose) * fvec`.
    qtf: Vec<T>,
    /// `m * nfree` column-major matrix holding the Jacobian (and its
    /// in-place QR factorization).
    fjac: JacMatrix<T>,
    /// Per-free-parameter lower/upper bounds (`None` means unbounded on that side).
    bounds: Vec<ParameterBounds<T>>,
    /// Cached: true if any free parameter has at least one limit.
    any_limit: bool,
    /// Model function to be evaluated.
    model: &'a mut U,
    /// Array of length n which defines the permutation matrix p such that `a * p = q * r`.
    /// Column j of p is column `ipvt(j)` of the identity matrix.
    ipvt: Vec<usize>,
    /// Per-parameter scale factors that form the diagonal matrix `D` in
    /// the LM damping term `(JᵀJ + λD²)`, length `npar`.
    diag: Vec<T>,
    /// Euclidean norm of the residual vector `fvec` at the current accepted point.
    fnorm: T,
    /// Euclidean norm of the residual vector `fvec` at the most recent trial point.
    /// Becomes `fnorm` once a step is accepted.
    fnorm1: T,
    /// Euclidean norm of the scaled parameter vector `diag * x`.
    xnorm: T,
    /// Variable which specifies an upper bound on the Euclidean norm of `d * x`.
    delta: T,
    /// Fit info.
    info: MPFitInfo,
    /// Squared residual norm at the starting parameters.
    orig_norm: T,
    /// Current Levenberg-Marquardt damping parameter (`(JᵀJ + λD²)`).
    lambda: T,
    /// Number of iterations of the algorithm.
    iter: usize,
    /// Fit configurations.
    cfg: &'a MPFitConfig<T>,
}

impl<'a, T, U> MPFit<'a, T, U>
where
    T: NumAssign + FloatConst,
    U: ImpedanceModel<T>,
{
    ///
    /// Build a fully-initialized `MPFit` ready to drive iterations.
    /// Validates the config, parses per-parameter bounds/fixed-vs-free
    /// status, and performs the initial model evaluation.
    ///
    pub fn try_new(
        model: &'a mut U,
        xall: Vec<T>,
        cfg: &'a MPFitConfig<T>,
    ) -> Result<Self, MPFitError> {
        let n_freqs = model.get_freqs().len();
        if n_freqs == 0 {
            return Err(MPFitError::Empty);
        }
        let m = 2 * n_freqs;
        let npar = xall.len();
        let mut fit = Self::new(m, npar, model, xall, cfg);
        fit.check_config()?;
        fit.parse_parameters()?;
        fit.init_lm()?;
        Ok(fit)
    }

    ///
    /// Build a new `MPFit` object.
    ///
    fn new(m: usize, npar: usize, model: &'a mut U, xall: Vec<T>, cfg: &'a MPFitConfig<T>) -> Self {
        Self {
            m,
            npar,
            nfree: 0,
            ifree: vec![],
            fvec: vec![T::zero(); m],
            nfev: 1,
            xnew: vec![T::zero(); npar],
            x: vec![],
            xall,
            qtf: vec![],
            fjac: JacMatrix::empty(),
            bounds: vec![],
            any_limit: false,
            model,
            ipvt: vec![0; npar],
            diag: vec![T::zero(); npar],
            fnorm: -T::one(),
            fnorm1: -T::one(),
            xnorm: -T::one(),
            delta: T::zero(),
            info: MPFitInfo::NotDone,
            orig_norm: T::zero(),
            lambda: T::zero(),
            iter: 1,
            cfg,
        }
    }

    ///
    /// Check the provided `MPFitConfig` for valid configurations.
    ///
    fn check_config(&self) -> Result<(), MPFitError> {
        if self.cfg.ftol <= T::zero()
            || self.cfg.xtol <= T::zero()
            || self.cfg.step_factor <= T::zero()
        {
            return Err(MPFitError::Input);
        }
        if self.m < self.nfree {
            return Err(MPFitError::DoF);
        }
        Ok(())
    }

    ///
    /// Parse and validate the model parameters.
    ///
    fn parse_parameters(&mut self) -> Result<(), MPFitError> {
        match &self.model.get_parameters() {
            None => {
                self.nfree = self.npar;
                self.ifree = (0..self.npar).collect();
                self.bounds = vec![ParameterBounds::default(); self.npar];
            }
            Some(pars) => {
                if pars.is_empty() {
                    return Err(MPFitError::Empty);
                }
                for (i, p) in pars.iter().enumerate() {
                    if !p.fit {
                        if p.bounds.lower.is_some_and(|x| self.xall[i] < x)
                            || p.bounds.upper.is_some_and(|x| self.xall[i] > x)
                        {
                            return Err(MPFitError::InitBounds);
                        }
                    } else {
                        if p.bounds.lower.is_some()
                            && p.bounds.upper.is_some()
                            && p.bounds.lower >= p.bounds.upper
                        {
                            return Err(MPFitError::Bounds);
                        }
                        self.nfree += 1;
                        self.ifree.push(i);
                        self.bounds.push(p.bounds);
                        if p.bounds.lower.is_some() || p.bounds.upper.is_some() {
                            self.any_limit = true;
                        }
                    }
                }
                if self.nfree == 0 {
                    return Err(MPFitError::NoFree);
                }
            }
        };
        if self.m < self.nfree {
            return Err(MPFitError::DoF);
        }
        Ok(())
    }

    ///
    /// Initialize Levenberg-Marquardt parameters and iteration counter.
    /// Starts by evaluating the model with the given starting parameters.
    ///
    fn init_lm(&mut self) -> Result<(), MPFitError> {
        self.model.evaluate(&self.xall, &mut self.fvec)?;
        self.nfev += 1;
        self.fnorm = self.fvec.enorm();
        self.orig_norm = self.fnorm * self.fnorm;
        self.xnew.copy_from_slice(&self.xall);
        self.x = self.ifree.iter().map(|&i| self.xall[i]).collect();
        self.qtf = vec![T::zero(); self.nfree];
        self.fjac = JacMatrix::zeros(self.m, self.nfree);

        Ok(())
    }

    ///
    /// Check if parameters are pegged at their upper/lower limits.
    ///
    pub fn check_limits(&mut self) {
        if self.any_limit {
            for j in 0..self.nfree {
                let bound = self.bounds[j];
                let lpegged = bound.lower.is_some_and(|l| self.x[j] == l);
                let upegged = bound.upper.is_some_and(|u| self.x[j] == u);
                let mut sum = T::zero();
                // If the parameter is pegged at a limit, compute the gradient direction
                if lpegged || upegged {
                    for (&f, &jc) in self.fvec.iter().zip(self.fjac.col(j)) {
                        sum += f * jc;
                    }
                }
                // If pegged at a limit and gradient is in the disallowed direction,
                // zero out the column.
                if (lpegged && sum > T::zero()) || (upegged && sum < T::zero()) {
                    self.fjac.col_mut(j).fill(T::zero());
                }
            }
        }
    }

    ///
    /// Calculate the Jacobian matrix.
    ///
    /// This function computes a forward-difference approximation
    /// to the m by n Jacobian matrix associated with a specified
    /// problem of m functions in n variables.
    ///
    pub fn fdjac2(&mut self) -> Result<(), MPFitError> {
        let eps = self.cfg.epsfcn.max(T::EPSILON).sqrt();
        self.fjac.fill(T::zero());
        let mut perturbed = vec![T::zero(); self.m];
        let mut ij = 0;
        // Any parameters requiring numerical derivatives
        for j in 0..self.nfree {
            let free_p = self.ifree[j];
            let temp = self.xnew[free_p];
            let mut h = eps * temp.abs();
            if h == T::zero() {
                h = eps;
            }
            if let Some(u) = self.bounds[j].upper {
                if temp > u - h {
                    h = -h;
                }
            }
            self.xnew[free_p] = temp + h;
            self.model.evaluate(&self.xnew, &mut perturbed)?;
            self.nfev += 1;
            self.xnew[free_p] = temp;
            for (&p, &fvec) in perturbed.iter().zip(&self.fvec) {
                self.fjac[ij] = (p - fvec) / h;
                ij += 1;
            }
        }
        Ok(())
    }

    ///
    /// Compute the QR factorization of the Jacobian.
    ///
    /// This function uses Householder transformations with column
    /// pivoting (optional) to compute a QR factorization of the
    /// m by n matrix a. That is, qrfac determines an orthogonal
    /// matrix q, a permutation matrix p, and an upper trapezoidal
    /// matrix r with diagonal elements of nonincreasing magnitude,
    /// such that a*p = q*r. the householder transformation for
    /// column k, k = 1,2,...,min(m,n), is of the form
    /// ```text
    ///     i - (1/u(k))*u*u
    /// ```
    /// where u has zeros in the first k-1 positions. The form of
    /// this transformation and the method of pivoting first
    /// appeared in the corresponding LINPACK subroutine.
    ///
    pub fn qrfac(&mut self, rdiag: &mut [T], acnorm: &mut [T]) {
        // `rdiag` receives the R diagonal (output, consumed by `transpose`).
        // `acnorm` receives the original column norms of A (output, consumed
        // by `scale`, `gnorm`, `rescale`).
        // `norms` is internal scratch used to detect when norms need
        // recomputing during pivoting.
        let mut norms = vec![T::zero(); self.nfree];
        // Compute the initial column norms and initialize several arrays.
        for (j, col) in self.fjac[0..self.m * self.nfree]
            .chunks_exact(self.m)
            .enumerate()
        {
            let norm = col.enorm();
            acnorm[j] = norm;
            rdiag[j] = norm;
            norms[j] = norm;
            self.ipvt[j] = j;
        }
        // Reduce a to r with householder transformations.
        for j in 0..self.m.min(self.nfree) {
            // Bring the column of largest norm into the pivot position.
            let mut kmax = j;
            for k in j + 1..self.nfree {
                if rdiag[k] > rdiag[kmax] {
                    kmax = k;
                }
            }
            if kmax != j {
                self.fjac.swap_cols(j, kmax);
                rdiag[kmax] = rdiag[j];
                norms[kmax] = norms[j];
                self.ipvt.swap(j, kmax);
            }
            let jj = j * (self.m + 1);
            let jjj = self.m - j + jj;
            let mut ajnorm = self.fjac[jj..jjj].enorm();
            if ajnorm == T::zero() {
                rdiag[j] = -ajnorm;
                continue;
            }
            if self.fjac[jj] < T::zero() {
                ajnorm = -ajnorm;
            }
            for fjac in self.fjac[jj..jjj].iter_mut() {
                *fjac /= ajnorm;
            }
            self.fjac[jj] += T::one();
            // Apply the transformation to the remaining columns
            // and update the norms.
            if j + 1 < self.nfree {
                for k in j + 1..self.nfree {
                    let mut sum = T::zero();
                    for i in j..self.m {
                        sum += self.fjac[self.m * j + i] * self.fjac[self.m * k + i];
                    }
                    let temp = sum / self.fjac[(j, j)];
                    for i in j..self.m {
                        let temp2 = self.fjac[self.m * j + i];
                        self.fjac[self.m * k + i] -= temp * temp2;
                    }
                    if rdiag[k] != T::zero() {
                        let temp = self.fjac[(j, k)] / rdiag[k];
                        let temp = (T::one() - temp.powi(2)).max(T::zero());
                        rdiag[k] *= temp.sqrt();
                        let temp = rdiag[k] / norms[k];
                        if temp.powi(2) * T::P05 < T::EPSILON {
                            let start = j + 1 + self.m * k;
                            rdiag[k] = self.fjac[start..start + self.m - j - 1].enorm();
                            norms[k] = rdiag[k];
                        }
                    }
                }
            }
            rdiag[j] = -ajnorm;
        }
    }

    ///
    /// On the first iteration and if user_scale is requested, scale according
    /// to the norms of the columns of the initial Jacobian,
    /// calculate the norm of the scaled x, and initialize the step bound delta.
    ///
    pub fn scale(&mut self, acnorm: &[T]) {
        if self.iter == 1 {
            if !self.cfg.do_user_scale {
                for j in 0..self.nfree {
                    self.diag[self.ifree[j]] = if acnorm[j] == T::zero() {
                        T::one()
                    } else {
                        acnorm[j]
                    };
                }
            }
            let mut scaled = vec![T::zero(); self.nfree];
            for j in 0..self.nfree {
                scaled[j] = self.diag[self.ifree[j]] * self.x[j];
            }
            self.xnorm = scaled.enorm();
            self.delta = self.cfg.step_factor * self.xnorm;
            if self.delta == T::zero() {
                self.delta = self.cfg.step_factor;
            }
        }
    }

    ///
    /// Fill xnew with x values for all free parameters.
    ///
    pub fn fill_xnew(&mut self) {
        for i in 0..self.nfree {
            self.xnew[self.ifree[i]] = self.x[i];
        }
    }

    ///
    /// Form (q transpose)*fvec and store the first n components in qtf.
    ///
    pub fn transpose(&mut self, rdiag: &[T]) {
        let mut qtfvec = self.fvec.clone();
        let mut jj = 0;
        for j in 0..self.nfree {
            let mut temp = self.fjac[jj];
            if temp != T::zero() {
                let mut sum = T::zero();
                let mut ij = jj;
                for i in j..self.m {
                    sum += self.fjac[ij] * qtfvec[i];
                    ij += 1;
                }
                temp = -sum / temp;
                ij = jj;
                for i in j..self.m {
                    qtfvec[i] += self.fjac[ij] * temp;
                    ij += 1;
                }
            }
            self.fjac[jj] = rdiag[j];
            jj += self.m + 1;
            self.qtf[j] = qtfvec[j];
        }
    }

    ///
    /// Check for overflow. This should be a cheap test here since FJAC
    /// has been reduced to a (small) square matrix, and the test is O(N^2).
    ///
    pub fn check_is_finite(&self) -> bool {
        self.cfg.no_finite_check || self.fjac.iter().all(|val| val.is_finite())
    }

    ///
    /// Compute the norm of the scaled gradient.
    ///
    pub fn gnorm(&self, acnorm: &[T]) -> T {
        let mut gnorm = T::zero();
        if self.fnorm != T::zero() {
            let mut jj = 0;
            for j in 0..self.nfree {
                let l = self.ipvt[j];
                if acnorm[l] != T::zero() {
                    let mut sum = T::zero();
                    for (ij, i) in (jj..).zip(0..=j) {
                        sum += self.fjac[ij] * (self.qtf[i] / self.fnorm);
                    }
                    gnorm = gnorm.max((sum / acnorm[l]).abs());
                }
                jj += self.m;
            }
        }
        gnorm
    }

    ///
    /// Terminate the fit and return the fit status.
    ///
    pub fn terminate(mut self) -> Result<MPFitStatus<T>, MPFitError> {
        for i in 0..self.nfree {
            self.xall[self.ifree[i]] = self.x[i];
        }
        // Compute number of pegged parameters
        let n_pegged = match self.model.get_parameters() {
            None => 0,
            Some(params) => {
                let mut n_pegged = 0;
                for (i, p) in params.iter().enumerate() {
                    if p.bounds.lower.is_some_and(|x| x == self.xall[i])
                        || p.bounds.upper.is_some_and(|x| x == self.xall[i])
                    {
                        n_pegged += 1;
                    }
                }
                n_pegged
            }
        };
        // Compute and return the covariance matrix and/or parameter errors
        self = self.covar();
        let mut covar = vec![T::zero(); self.npar * self.npar];
        for j in 0..self.nfree {
            let k = self.ifree[j] * self.npar;
            for i in 0..self.nfree {
                covar[k + self.ifree[i]] = self.fjac[(i, j)];
            }
        }
        let mut xerror = vec![T::zero(); self.npar];
        for j in 0..self.nfree {
            let cc = self.fjac[(j, j)];
            if cc > T::zero() {
                xerror[self.ifree[j]] = cc.sqrt();
            }
        }
        let best_norm = self.fnorm.max(self.fnorm1);
        Ok(MPFitStatus {
            info: self.info,
            best_norm: best_norm * best_norm,
            orig_norm: self.orig_norm,
            n_iter: self.iter,
            n_fev: self.nfev,
            n_par: self.npar,
            n_free: self.nfree,
            n_pegged,
            n_func: self.m,
            residuals: self.fvec,
            xerror,
            covar,
            x: self.xall,
        })
    }

    ///
    /// Calculate the covariance matrix.
    ///
    /// Given an m by n matrix a, the problem is to determine
    /// the covariance matrix corresponding to a, defined as
    /// ```text
    ///       inverse(a *a)
    /// ```
    /// This subroutine completes the solution of the problem
    /// if it is provided with the necessary information from the
    /// QR factorization, with column pivoting, of a. That is, if
    /// a*p = q*r, where p is a permutation matrix, q has orthogonal
    /// columns, and r is an upper triangular matrix with diagonal
    /// elements of nonincreasing magnitude, then covar expects
    /// the full upper triangle of r and the permutation matrix p.
    /// The covariance matrix is then computed as
    /// ```text
    ///       p*inverse(r *r)*p
    /// ```
    /// If a is nearly rank deficient, it may be desirable to compute
    /// the covariance matrix corresponding to the linearly independent
    /// columns of a. To define the numerical rank of a, covar uses
    /// the tolerance tol. If l is the largest integer such that
    /// ```text
    ///       abs(r(l,l)) > tol*abs(r(1,1))
    /// ```
    /// then covar computes the covariance matrix corresponding to
    /// the first l columns of r. For k greater than l, column
    /// and row ipvt(k) of the covariance matrix are set to zero.
    ///
    pub fn covar(mut self) -> Self {
        // Form the inverse of r in the full upper triangle of r.
        let tolr = self.cfg.covtol * self.fjac[(0, 0)].abs();
        let mut l: isize = -1;
        for k in 0..self.nfree {
            if self.fjac[(k, k)].abs() <= tolr {
                break;
            }
            self.fjac[(k, k)] = T::one() / self.fjac[(k, k)];
            for j in 0..k {
                let temp = self.fjac[(k, k)] * self.fjac[(j, k)];
                self.fjac[(j, k)] = T::zero();
                for i in 0..=j {
                    let dfjac = -temp * self.fjac[(i, j)];
                    self.fjac[(i, k)] += dfjac;
                }
            }
            l = k as isize;
        }
        // Form the full upper triangle of the inverse of (r transpose)*r
        // in the full upper triangle of r
        if l >= 0 {
            let l = l as usize;
            for k in 0..=l {
                for j in 0..k {
                    let temp = self.fjac[(j, k)];
                    for i in 0..=j {
                        let dfjac = temp * self.fjac[(i, k)];
                        self.fjac[(i, j)] += dfjac;
                    }
                }
                let temp = self.fjac[(k, k)];
                for i in 0..=k {
                    self.fjac[(i, k)] *= temp;
                }
            }
        }
        // For the full lower triangle of the covariance matrix
        // in the strict lower triangle of r
        let mut covar_diag = vec![T::zero(); self.nfree];
        for j in 0..self.nfree {
            let jj = self.ipvt[j];
            let sing = j as isize > l;
            for i in 0..=j {
                if sing {
                    self.fjac[(i, j)] = T::zero();
                }
                let ii = self.ipvt[i];
                if ii > jj {
                    self.fjac[(ii, jj)] = self.fjac[(i, j)];
                }
                if ii < jj {
                    self.fjac[(jj, ii)] = self.fjac[(i, j)];
                }
            }
            covar_diag[jj] = self.fjac[(j, j)];
        }
        // Symmetrize the covariance matrix in r
        for j in 0..self.nfree {
            for i in 0..j {
                self.fjac[(i, j)] = self.fjac[(j, i)];
            }
            self.fjac[(j, j)] = covar_diag[j];
        }
        self
    }

    pub fn rescale(&mut self, acnorm: &[T]) {
        if !self.cfg.do_user_scale {
            for j in 0..self.nfree {
                let i = self.ifree[j];
                self.diag[i] = self.diag[i].max(acnorm[j]);
            }
        }
    }

    ///
    /// Given an m by nfree matrix a, an nfree by nfree nonsingular diagonal
    /// matrix d, an m-vector b, and a positive number delta,
    /// the problem is to determine a value for the parameter
    /// par such that if `step` solves the system
    /// ```text
    ///     a*step = b ,   sqrt(par)*d*step = 0
    /// ```
    /// in the least squares sense, and dxnorm is the Euclidean
    /// norm of d*step, then either par is zero and
    /// ```text
    ///     (dxnorm-delta) < 0.1*delta
    /// ```
    /// or par is positive and
    /// ```text
    ///     abs(dxnorm-delta) < 0.1*delta
    /// ```
    /// This subroutine completes the solution of the problem
    /// if it is provided with the necessary information from the
    /// QR factorization, with column pivoting, of a. That is, if
    /// a*p = q*fjack, where p is a permutation matrix, q has orthogonal
    /// columns, and fjack is an upper triangular matrix with diagonal
    /// elements of nonincreasing magnitude, then lmpar expects
    /// the full upper triangle of fjack, the permutation matrix p,
    /// and the first nfree components of (q transpose)*b. On output
    /// lmpar also provides an upper triangular matrix s such that
    /// ```text
    ///     p *(a *a + par*d*d)*p = s *s
    /// ```
    /// s is employed within lmpar and may be of separate interest.
    ///
    /// Only a few iterations are generally needed for convergence
    /// of the algorithm. If, however, the limit of 10 iterations
    /// is reached, then the output par will contain the best
    /// value obtained so far.
    ///
    pub fn lmpar(&mut self, step: &mut [T]) {
        let mut work = vec![T::zero(); self.nfree];
        let mut scaled_step = vec![T::zero(); self.nfree];
        let mut sdiag = vec![T::zero(); self.nfree];

        let nsing = self.gauss_newton_step(step, &mut work);

        for j in 0..self.nfree {
            scaled_step[j] = self.diag[self.ifree[j]] * step[j];
        }
        let mut dxnorm = scaled_step.enorm();
        let mut fp = dxnorm - self.delta;
        if fp <= self.delta * T::P1 {
            self.lambda = T::zero();
            return;
        }

        let mut lambda_lower =
            self.compute_lambda_lower(nsing, dxnorm, &scaled_step, fp, &mut work);
        let (gnorm, mut lambda_upper) = self.compute_lambda_upper(&mut work);

        // Clamp the input par into the bracket. Note: the original CMPFIT
        // port uses `.max(lambda_upper)` on both sides — preserved here to keep
        // bit-identical behavior with the reference implementation.
        self.lambda = self.lambda.max(lambda_lower);
        self.lambda = self.lambda.max(lambda_upper);
        if self.lambda == T::zero() {
            self.lambda = gnorm / dxnorm;
        }

        for iter in 1..=10 {
            if self.lambda == T::zero() {
                self.lambda = T::MIN_POSITIVE.max(lambda_upper * T::P0001);
            }
            let temp = self.lambda.sqrt();
            for j in 0..self.nfree {
                work[j] = temp * self.diag[self.ifree[j]];
            }
            self.qrsolv(&work, step, &mut sdiag);
            for j in 0..self.nfree {
                scaled_step[j] = self.diag[self.ifree[j]] * step[j];
            }
            dxnorm = scaled_step.enorm();
            let prev_fp = fp;
            fp = dxnorm - self.delta;
            if fp.abs() <= self.delta * T::P1
                || (lambda_lower == T::zero() && fp <= prev_fp && prev_fp < T::zero())
                || iter >= 10
            {
                return;
            }
            let parc = self.par_correction(dxnorm, fp, &scaled_step, &sdiag, &mut work);
            if fp > T::zero() {
                lambda_lower = lambda_lower.max(self.lambda);
            }
            if fp < T::zero() {
                lambda_upper = lambda_upper.min(self.lambda);
            }
            self.lambda = lambda_lower.max(self.lambda + parc);
        }
    }

    ///
    /// Compute the (least-squares) Gauss-Newton direction by
    /// back-substitution of the QR factorization. Stores the permuted
    /// result in `step`; uses `work` as length-`nfree` scratch. Returns
    /// the smallest column index whose R diagonal is zero, or `nfree` if
    /// R is non-singular.
    ///
    fn gauss_newton_step(&self, step: &mut [T], work: &mut [T]) -> usize {
        let mut nsing = self.nfree;
        for j in 0..self.nfree {
            work[j] = self.qtf[j];
            if self.fjac[(j, j)] == T::zero() && nsing == self.nfree {
                nsing = j;
            }
            if nsing < self.nfree {
                work[j] = T::zero();
            }
        }
        for j in (0..nsing).rev() {
            work[j] /= self.fjac[(j, j)];
            let temp = work[j];
            for i in 0..j {
                work[i] -= self.fjac[self.m * j + i] * temp;
            }
        }
        for j in 0..self.nfree {
            step[self.ipvt[j]] = work[j];
        }
        nsing
    }

    ///
    /// Compute the lower bracket `lambda_lower` on the LM parameter. Returns zero
    /// when the Jacobian is rank-deficient (no useful lower bound).
    ///
    fn compute_lambda_lower(
        &self,
        nsing: usize,
        dxnorm: T,
        scaled_step: &[T],
        fp: T,
        work: &mut [T],
    ) -> T {
        if nsing < self.nfree {
            return T::zero();
        }
        self.newton_correction(dxnorm, scaled_step, work);
        let mut jj = 0;
        for j in 0..self.nfree {
            let mut sum = T::zero();
            if j > 0 {
                let mut ij = jj;
                for i in 0..j {
                    sum += self.fjac[ij] * work[i];
                    ij += 1;
                }
            }
            work[j] = (work[j] - sum) / self.fjac[(j, j)];
            jj += self.m;
        }
        let temp = work.enorm();
        ((fp / self.delta) / temp) / temp
    }

    ///
    /// Compute the upper bracket `lambda_upper` on the LM parameter from the
    /// scaled gradient norm. Returns `(gnorm, lambda_upper)`.
    ///
    fn compute_lambda_upper(&self, work: &mut [T]) -> (T, T) {
        let mut jj = 0;
        for j in 0..self.nfree {
            let mut sum = T::zero();
            let mut ij = jj;
            for i in 0..=j {
                sum += self.fjac[ij] * self.qtf[i];
                ij += 1;
            }
            let l = self.ipvt[j];
            work[j] = sum / self.diag[self.ifree[l]];
            jj += self.m;
        }
        let gnorm = work.enorm();
        let mut lambda_upper = gnorm / self.delta;
        if lambda_upper == T::zero() {
            lambda_upper = T::MIN_POSITIVE / self.delta.min(T::P1);
        }
        (gnorm, lambda_upper)
    }

    ///
    /// Compute the iterative correction `parc` to the LM parameter using
    /// the s-diagonal returned by qrsolv. Modifies `work` in place.
    ///
    fn par_correction(
        &self,
        dxnorm: T,
        fp: T,
        scaled_step: &[T],
        sdiag: &[T],
        work: &mut [T],
    ) -> T {
        self.newton_correction(dxnorm, scaled_step, work);
        let mut jj = 0;
        for j in 0..self.nfree {
            work[j] /= sdiag[j];
            let temp = work[j];
            if j + 1 < self.nfree {
                let mut ij = j + 1 + jj;
                for i in j + 1..self.nfree {
                    work[i] -= self.fjac[ij] * temp;
                    ij += 1;
                }
            }
            jj += self.m;
        }
        let temp = work.enorm();
        ((fp / self.delta) / temp) / temp
    }

    ///
    /// Compute the Newton correction.
    ///
    fn newton_correction(&self, dxnorm: T, scaled_step: &[T], work: &mut [T]) {
        for j in 0..self.nfree {
            let l = self.ipvt[j];
            work[j] = self.diag[self.ifree[l]] * (scaled_step[l] / dxnorm);
        }
    }

    ///
    /// Given an m by n matrix a, an n by n diagonal matrix d,
    /// and an m-vector b, the problem is to determine an x which
    /// solves the system
    /// ```text
    /// a*x = b ,     d*x = 0
    /// ```
    /// in the least squares sense.
    ///
    /// This subroutine completes the solution of the problem
    /// if it is provided with the necessary information from the
    /// QR factorization, with column pivoting, of a. That is, if
    /// a*p = q*r, where p is a permutation matrix, q has orthogonal
    /// columns, and r is an upper triangular matrix with diagonal
    /// elements of nonincreasing magnitude, then qrsolv expects
    /// the full upper triangle of r, the permutation matrix p,
    /// and the first n components of (q transpose)*b. The system
    /// a*x = b, d*x = 0, is then equivalent to
    /// ```text
    /// r*z = q *b ,  p *d*p*z = 0 ,
    /// ```
    /// where x = p*z. If this system does not have full rank,
    /// then a least squares solution is obtained. On output qrsolv
    /// also provides an upper triangular matrix s such that
    /// ```text
    /// p *(a *a + d*d)*p = s *s .
    /// ```
    /// s is computed within qrsolv and may be of separate interest.
    ///
    fn qrsolv(&mut self, dpar: &[T], step: &mut [T], sdiag: &mut [T]) {
        // `dpar` is the input sqrt(par)*diag (length nfree).
        // `step` is the output LM step (length nfree), permuted back at end.
        // `sdiag` is the output s-diagonal (length nfree), used by lmpar.
        // `qtb` (length nfree) is local: a working copy of qtf that becomes z.
        // `r_diag` (length nfree) is local: temporary save of R diagonal which
        // is restored at the end so `fjac` is unchanged on exit.
        let mut qtb = vec![T::zero(); self.nfree];
        let mut r_diag = vec![T::zero(); self.nfree];
        // Copy r and (q transpose)*b to preserve input and initialize s.
        // in particular, save the diagonal elements of r in r_diag.
        let mut kk = 0;
        for j in 0..self.nfree {
            // Mirror row j of r into column j
            for i in j..self.nfree {
                self.fjac[self.m * j + i] = self.fjac[j + self.m * i];
            }
            r_diag[j] = self.fjac[kk];
            qtb[j] = self.qtf[j];
            kk += self.m + 1;
        }
        // Eliminate the diagonal matrix d using a givens rotation.
        for j in 0..self.nfree {
            // Prepare the row of d to be eliminated, locating the
            // diagonal element using p from the qr factorization.
            let l = self.ipvt[j];
            if dpar[l] != T::zero() {
                for k in j..self.nfree {
                    sdiag[k] = T::zero();
                }
                sdiag[j] = dpar[l];
                // The transformations to eliminate the row of d
                // modify only a single element of (q transpose)*b
                // beyond the first n, which is initially zero.
                let mut qtbpj = T::zero();
                for k in j..self.nfree {
                    // Determine a givens rotation which eliminates the
                    // appropriate element in the current row of d.
                    if sdiag[k] == T::zero() {
                        continue;
                    }
                    let kk = k + self.m * k;
                    let (sinx, cosx) = if self.fjac[kk].abs() < sdiag[k].abs() {
                        let cotan = self.fjac[kk] / sdiag[k];
                        let sinx = T::HALF / (cotan * cotan * T::P25 + T::P25).sqrt();
                        let cosx = sinx * cotan;
                        (sinx, cosx)
                    } else {
                        let tanx = sdiag[k] / self.fjac[kk];
                        let cosx = T::HALF / (tanx * tanx * T::P25 + T::P25).sqrt();
                        let sinx = cosx * tanx;
                        (sinx, cosx)
                    };
                    // Compute the modified diagonal element of r and
                    // the modified element of ((q transpose)*b,0).
                    self.fjac[kk] = cosx * self.fjac[kk] + sinx * sdiag[k];
                    let temp = cosx * qtb[k] + sinx * qtbpj;
                    qtbpj = -sinx * qtb[k] + cosx * qtbpj;
                    qtb[k] = temp;
                    // Accumulate the transformation in the row of s.
                    for i in k + 1..self.nfree {
                        let j = self.m * k + i;
                        let f = self.fjac[j];
                        let s = sdiag[i];
                        self.fjac[j] = cosx * f + sinx * s;
                        sdiag[i] = cosx * s - sinx * f;
                    }
                }
            }
            // Store the diagonal element of s and restore
            // the corresponding diagonal element of r.
            let kk = j + self.m * j;
            sdiag[j] = self.fjac[kk];
            self.fjac[kk] = r_diag[j];
        }
        // Solve the triangular system for z. if the system is
        // singular, then obtain a least squares solution.
        let mut nsing = self.nfree;
        for j in 0..self.nfree {
            if sdiag[j] == T::zero() && nsing == self.nfree {
                nsing = j;
            }
            if nsing < self.nfree {
                qtb[j] = T::zero();
            }
        }
        for j in (0..nsing).rev() {
            let mut sum = T::zero();
            for i in j + 1..nsing {
                sum += self.fjac[self.m * j + i] * qtb[i];
            }
            qtb[j] = (qtb[j] - sum) / sdiag[j];
        }
        // Permute the components of z back to components of x.
        for j in 0..self.nfree {
            step[self.ipvt[j]] = qtb[j];
        }
    }

    ///
    /// Take one Levenberg-Marquardt trial step. Returns `true` if the step was
    /// accepted (caller should re-evaluate the Jacobian), `false` if it was
    /// rejected (caller should shrink the trust region and try again). Either
    /// way, the caller must also inspect `self.info` — if it is not
    /// `NotDone`, the fit has terminated and `terminate()` should be called.
    ///
    pub fn iterate(&mut self, gnorm: T, step: &mut [T]) -> Result<bool, MPFitError> {
        let mut trial_x = vec![T::zero(); self.nfree];
        let mut work = vec![T::zero(); self.nfree];
        let mut trial_residuals = vec![T::zero(); self.m];

        for j in 0..self.nfree {
            step[j] = -step[j];
        }

        let alpha = self.apply_bounds(step, &mut trial_x);

        for j in 0..self.nfree {
            work[j] = self.diag[self.ifree[j]] * step[j];
        }
        let pnorm = work.enorm();
        if self.iter == 1 {
            self.delta = self.delta.min(pnorm);
        }

        for i in 0..self.nfree {
            self.xnew[self.ifree[i]] = trial_x[i];
        }
        self.fnorm1 = self.evaluate_trial(&mut trial_residuals)?;

        let (actred, prered, dirder, ratio) =
            self.compute_reductions(step, alpha, pnorm, &mut work);

        self.update_trust_region(actred, dirder, ratio, pnorm);

        let accepted = ratio >= T::P0001;
        if accepted {
            for j in 0..self.nfree {
                self.x[j] = trial_x[j];
                work[j] = self.diag[self.ifree[j]] * self.x[j];
            }
            self.fvec.copy_from_slice(&trial_residuals);
            self.xnorm = work.enorm();
            self.fnorm = self.fnorm1;
            self.iter += 1;
        }

        self.check_convergence(actred, prered, ratio, gnorm);

        Ok(accepted)
    }

    ///
    /// Apply parameter limits to `step` and produce the candidate next
    /// position `trial_x`. Returns `alpha`, the fraction of the full step
    /// actually taken (`1.0` when no bound is hit).
    ///
    fn apply_bounds(&self, step: &mut [T], trial_x: &mut [T]) -> T {
        let mut alpha = T::one();
        if !self.any_limit {
            for j in 0..self.nfree {
                trial_x[j] = self.x[j] + step[j];
            }
            return alpha;
        }
        // Respect the limits. If the step would go out of bounds, shorten
        // it so it lands exactly on the limit.
        for j in 0..self.nfree {
            let bound = self.bounds[j];
            let lpegged = bound.lower.is_some_and(|l| self.x[j] <= l);
            let upegged = bound.upper.is_some_and(|u| self.x[j] >= u);
            let dstep = step[j].abs() > T::EPSILON;
            if lpegged && step[j] < T::zero() {
                step[j] = T::zero();
            }
            if upegged && step[j] > T::zero() {
                step[j] = T::zero();
            }
            if dstep {
                if let Some(l) = bound.lower {
                    if self.x[j] + step[j] < l {
                        alpha = alpha.min((l - self.x[j]) / step[j]);
                    }
                }
                if let Some(u) = bound.upper {
                    if self.x[j] + step[j] > u {
                        alpha = alpha.min((u - self.x[j]) / step[j]);
                    }
                }
            }
        }
        for j in 0..self.nfree {
            step[j] *= alpha;
            trial_x[j] = self.x[j] + step[j];
            // Snap to the limit when within EPSILON of it, so a step that
            // intended to land on the boundary does so exactly.
            let bound = self.bounds[j];
            if let Some(u) = bound.upper {
                let sgnu = if u >= T::zero() { T::one() } else { -T::one() };
                let ulim1 = u * (T::one() - sgnu * T::EPSILON)
                    - if u == T::zero() {
                        T::EPSILON
                    } else {
                        T::zero()
                    };
                if trial_x[j] >= ulim1 {
                    trial_x[j] = u;
                }
            }
            if let Some(l) = bound.lower {
                let sgnl = if l >= T::zero() { T::one() } else { -T::one() };
                let llim1 = l * (T::one() + sgnl * T::EPSILON)
                    + if l == T::zero() {
                        T::EPSILON
                    } else {
                        T::zero()
                    };
                if trial_x[j] <= llim1 {
                    trial_x[j] = l;
                }
            }
        }
        alpha
    }

    ///
    /// Evaluate the model at `self.xnew`, store the residual in
    /// `trial_residualsuals`, and return its Euclidean norm.
    ///
    fn evaluate_trial(&mut self, trial_residualsuals: &mut [T]) -> Result<T, MPFitError> {
        self.model.evaluate(&self.xnew, trial_residualsuals)?;
        self.nfev += 1;
        Ok(trial_residualsuals.enorm())
    }

    ///
    /// Compute the scaled actual / predicted reductions and directional
    /// derivative. Returns `(actred, prered, dirder, ratio)`. `alpha` is
    /// the fraction of the full LM step actually taken (set by
    /// `apply_bounds`).
    ///
    fn compute_reductions(&self, step: &[T], alpha: T, pnorm: T, work: &mut [T]) -> (T, T, T, T) {
        let actred = if self.fnorm1 * T::P1 < self.fnorm {
            let temp = self.fnorm1 / self.fnorm;
            T::one() - temp * temp
        } else {
            -T::one()
        };
        let mut jj = 0;
        for j in 0..self.nfree {
            work[j] = T::zero();
            let l = self.ipvt[j];
            let temp = step[l];
            let mut ij = jj;
            for i in 0..=j {
                work[i] += self.fjac[ij] * temp;
                ij += 1;
            }
            jj += self.m;
        }
        let temp1 = work.enorm() * alpha / self.fnorm;
        let temp2 = ((alpha * self.lambda).sqrt() * pnorm) / self.fnorm;
        let temp11 = temp1 * temp1;
        let temp22 = temp2 * temp2;
        let prered = temp11 + temp22 / T::HALF;
        let dirder = -(temp11 + temp22);
        let ratio = if prered != T::zero() {
            actred / prered
        } else {
            T::zero()
        };
        (actred, prered, dirder, ratio)
    }

    ///
    /// Update the trust-region radius `delta` and LM damping `lambda` from
    /// the reduction ratio.
    ///
    fn update_trust_region(&mut self, actred: T, dirder: T, ratio: T, pnorm: T) {
        if ratio <= T::P25 {
            let mut temp: T = if actred >= T::zero() {
                T::HALF
            } else {
                dirder * T::HALF / (dirder + actred * T::HALF)
            };
            if self.fnorm1 * T::P1 >= self.fnorm || temp < T::P1 {
                temp = T::P1;
            }
            self.delta = temp * self.delta.min(pnorm / T::P1);
            self.lambda /= temp;
        } else if self.lambda == T::zero() || ratio >= T::P75 {
            self.delta = pnorm / T::HALF;
            self.lambda *= T::HALF;
        }
    }

    ///
    /// Set `self.info` based on convergence and termination conditions.
    /// The stringent-tolerance block only runs when no Convergence* status
    /// has already been set, so it cannot downgrade a converged fit to
    /// MaxIterReached or *NoImprovement.
    ///
    fn check_convergence(&mut self, actred: T, prered: T, ratio: T, gnorm: T) {
        if actred.abs() <= self.cfg.ftol && prered <= self.cfg.ftol && ratio * T::HALF <= T::one() {
            self.info = MPFitInfo::ConvergenceChi;
        }
        if self.delta <= self.cfg.xtol * self.xnorm {
            self.info = MPFitInfo::ConvergencePar;
        }
        if actred.abs() <= self.cfg.ftol
            && prered <= self.cfg.ftol
            && ratio * T::HALF <= T::one()
            && self.info == MPFitInfo::ConvergencePar
        {
            self.info = MPFitInfo::ConvergenceBoth;
        }
        if self.info == MPFitInfo::NotDone {
            if self.cfg.max_fev > 0 && self.nfev >= self.cfg.max_fev {
                self.info = MPFitInfo::MaxIterReached;
            }
            if self.iter >= self.cfg.max_iter {
                self.info = MPFitInfo::MaxIterReached;
            }
            if actred.abs() <= T::EPSILON && prered <= T::EPSILON && ratio * T::HALF <= T::one() {
                self.info = MPFitInfo::FtolNoImprovement;
            }
            if self.delta <= T::EPSILON * self.xnorm {
                self.info = MPFitInfo::XtolNoImprovement;
            }
            if gnorm <= T::EPSILON {
                self.info = MPFitInfo::GtolNoImprovement;
            }
        }
    }

    ///
    /// Check for convergence in orthogonality.
    ///
    pub fn check_convergence_ortho(&mut self, acnorm: &[T]) {
        if self.gnorm(acnorm) <= self.cfg.gtol {
            self.info = MPFitInfo::ConvergenceDir;
        }
    }

    ///
    /// Check if `max_iter == 0` and set `self.info` if so.
    ///
    pub fn check_no_iter(&mut self) {
        if self.cfg.max_iter == 0 {
            self.info = MPFitInfo::MaxIterReached;
        }
    }

    ///
    /// True if `self.info` is any status other than not done.
    ///
    pub fn is_done(&self) -> bool {
        self.info != MPFitInfo::NotDone
    }

    ///
    /// Number of free parameters.
    ///
    pub fn nfree(&self) -> usize {
        self.nfree
    }
}
