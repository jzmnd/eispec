//!
//! Copied and modified from rmpfit crate Copyright (c) Vadim Dyadkin
//! Rust implementation of [CMPFIT](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html)
//!
use crate::fit::config::MPFitConfig;
use crate::fit::status::MPFitStatus;
use num::traits::NumAssign;

use crate::constants::FloatConst;
use crate::fit::enorm::ENorm;
use crate::fit::enums::{MPFitDone, MPFitError, MPFitInfo};
use crate::fit::ImpedanceDataFitter;

pub struct MPFit<'a, T, F> {
    pub m: usize,
    pub npar: usize,
    pub nfree: usize,
    pub ifree: Vec<usize>,
    pub fvec: Vec<T>,
    pub nfev: usize,
    pub xnew: Vec<T>,
    pub x: Vec<T>,
    pub xall: &'a mut [T],
    pub qtf: Vec<T>,
    pub fjac: Vec<T>,
    pub step: Vec<T>,
    pub dstep: Vec<T>,
    pub qllim: Vec<bool>,
    pub qulim: Vec<bool>,
    pub llim: Vec<T>,
    pub ulim: Vec<T>,
    pub qanylim: bool,
    pub f: &'a mut F,
    pub wa1: Vec<T>,
    pub wa2: Vec<T>,
    pub wa3: Vec<T>,
    pub wa4: Vec<T>,
    pub ipvt: Vec<usize>,
    pub diag: Vec<T>,
    pub fnorm: T,
    pub fnorm1: T,
    pub xnorm: T,
    pub delta: T,
    pub info: MPFitInfo,
    pub orig_norm: T,
    pub par: T,
    pub iter: usize,
    pub cfg: &'a MPFitConfig<T>,
}

impl<'a, T, F> MPFit<'a, T, F>
where
    T: NumAssign + FloatConst,
    F: ImpedanceDataFitter<T>,
{
    pub fn new(
        f: &'a mut F,
        xall: &'a mut [T],
        cfg: &'a MPFitConfig<T>,
    ) -> Result<Self, MPFitError> {
        let m = f.freqs().len();
        let npar = xall.len();
        if m == 0 {
            Err(MPFitError::Empty)
        } else {
            Ok(Self {
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
                fjac: vec![],
                step: vec![],
                dstep: vec![],
                qllim: vec![],
                qulim: vec![],
                llim: vec![],
                ulim: vec![],
                qanylim: false,
                f,
                wa1: vec![T::zero(); npar],
                wa2: vec![T::zero(); m],
                wa3: vec![T::zero(); npar],
                wa4: vec![T::zero(); m],
                ipvt: vec![0; npar],
                diag: vec![T::zero(); npar],
                fnorm: -T::one(),
                fnorm1: -T::one(),
                xnorm: -T::one(),
                delta: T::zero(),
                info: MPFitInfo::NotDone,
                orig_norm: T::zero(),
                par: T::zero(),
                iter: 1,
                cfg,
            })
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
        let mut ij = 0;
        // Any parameters requiring numerical derivatives
        for j in 0..self.nfree {
            let free_p = self.ifree[j];
            let temp = self.xnew[free_p];
            let mut h = eps * temp.abs();
            if free_p < self.step.len() && self.step[free_p] > T::zero() {
                h = self.step[free_p];
            }
            if free_p < self.dstep.len() && self.dstep[free_p] > T::zero() {
                h = (self.dstep[free_p] * temp).abs();
            }
            if h == T::zero() {
                h = eps;
            }
            if j < self.qulim.len()
                && self.qulim[j]
                && j < self.ulim.len()
                && temp > self.ulim[j] - h
            {
                h = -h;
            }
            self.xnew[free_p] = temp + h;
            self.f.evaluate(&self.xnew, &mut self.wa4)?;
            self.nfev += 1;
            self.xnew[free_p] = temp;
            for (wa4, fvec) in self.wa4.iter().zip(&self.fvec) {
                self.fjac[ij] = (*wa4 - *fvec) / h;
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
    pub fn qrfac(&mut self) {
        // Compute the initial column norms and initialize several arrays.
        for (j, ij) in (0..self.nfree).zip((0..self.m * self.nfree).step_by(self.m)) {
            self.wa2[j] = self.fjac[ij..ij + self.m].enorm();
            self.wa1[j] = self.wa2[j];
            self.wa3[j] = self.wa1[j];
            self.ipvt[j] = j;
        }
        // Reduce a to r with householder transformations.
        for j in 0..self.m.min(self.nfree) {
            // Bring the column of largest norm into the pivot position.
            let mut kmax = j;
            for k in j..self.nfree {
                if self.wa1[k] > self.wa1[kmax] {
                    kmax = k;
                }
            }
            if kmax != j {
                let mut ij = self.m * j;
                let mut jj = self.m * kmax;
                for _ in 0..self.m {
                    self.fjac.swap(jj, ij);
                    ij += 1;
                    jj += 1;
                }
                self.wa1[kmax] = self.wa1[j];
                self.wa3[kmax] = self.wa3[j];
                self.ipvt.swap(j, kmax);
            }
            let jj = j + self.m * j;
            let jjj = self.m - j + jj;
            let mut ajnorm = self.fjac[jj..jjj].enorm();
            if ajnorm == T::zero() {
                self.wa1[j] = -ajnorm;
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
            let jp1 = j + 1;
            if jp1 < self.nfree {
                for k in jp1..self.nfree {
                    let mut sum = T::zero();
                    let mut ij = j + self.m * k;
                    let mut jj = j + self.m * j;
                    for _ in j..self.m {
                        sum += self.fjac[jj] * self.fjac[ij];
                        ij += 1;
                        jj += 1;
                    }
                    let temp = sum / self.fjac[j + self.m * j];
                    ij = j + self.m * k;
                    jj = j + self.m * j;
                    for _ in j..self.m {
                        let dfjac = temp * self.fjac[jj];
                        self.fjac[ij] -= dfjac;
                        ij += 1;
                        jj += 1;
                    }
                    if self.wa1[k] != T::zero() {
                        let temp = self.fjac[j + self.m * k] / self.wa1[k];
                        let temp = (T::one() - temp.powi(2)).max(T::zero());
                        self.wa1[k] *= temp.sqrt();
                        let temp = self.wa1[k] / self.wa3[k];
                        if temp.powi(2) * T::P05 < T::EPSILON {
                            let start = jp1 + self.m * k;
                            self.wa1[k] = self.fjac[start..start + self.m - j - 1].enorm();
                            self.wa3[k] = self.wa1[k];
                        }
                    }
                }
            }
            self.wa1[j] = -ajnorm;
        }
    }

    ///
    /// Parse and validate the model parameters.
    ///
    pub fn parse_params(&mut self) -> Result<(), MPFitError> {
        match &self.f.model_params() {
            None => {
                self.nfree = self.npar;
                self.ifree = (0..self.npar).collect();
            }
            Some(pars) => {
                if pars.is_empty() {
                    return Err(MPFitError::Empty);
                }
                for (i, p) in pars.iter().enumerate() {
                    if !p.fit {
                        if p.limit_lower.is_some_and(|x| self.xall[i] < x)
                            || p.limit_upper.is_some_and(|x| self.xall[i] > x)
                        {
                            return Err(MPFitError::InitBounds);
                        }
                    } else {
                        if p.limit_lower.is_some()
                            && p.limit_upper.is_some()
                            && p.limit_lower >= p.limit_upper
                        {
                            return Err(MPFitError::Bounds);
                        }
                        self.nfree += 1;
                        self.ifree.push(i);
                        self.qllim.push(p.limit_lower.is_some());
                        self.qulim.push(p.limit_upper.is_some());
                        self.llim.push(p.limit_lower.unwrap_or(T::zero()));
                        self.ulim.push(p.limit_upper.unwrap_or(T::zero()));
                        if p.limit_lower.is_some() || p.limit_upper.is_some() {
                            self.qanylim = true;
                        }
                    }
                    self.step.push(T::zero());
                    self.dstep.push(T::zero());
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
    ///
    pub fn init_lm(&mut self) -> Result<(), MPFitError> {
        self.f.evaluate(self.xall, &mut self.fvec)?;
        self.nfev += 1;
        self.fnorm = self.fvec.enorm();
        self.orig_norm = self.fnorm * self.fnorm;
        self.xnew.copy_from_slice(self.xall);
        self.x = self.ifree.iter().map(|&i| self.xall[i]).collect();
        self.qtf = vec![T::zero(); self.nfree];
        self.fjac = vec![T::zero(); self.m * self.nfree];

        Ok(())
    }

    ///
    /// Check if parameters are pegged at their upper/lower limits.
    ///
    pub fn check_limits(&mut self) {
        if self.qanylim {
            for j in 0..self.nfree {
                let lpegged = j < self.qllim.len() && self.x[j] == self.llim[j];
                let upegged = j < self.qulim.len() && self.x[j] == self.ulim[j];
                let mut sum = T::zero();
                // If the parameter is pegged at a limit, compute the gradient direction
                let ij = j * self.m;
                if lpegged || upegged {
                    for i in 0..self.m {
                        sum += self.fvec[i] * self.fjac[ij + i];
                    }
                }
                // If pegged at lower limit and gradient is toward negative then
                // reset gradient to zero
                if lpegged && sum > T::zero() {
                    for i in 0..self.m {
                        self.fjac[ij + i] = T::zero();
                    }
                }
                // If pegged at upper limit and gradient is toward positive then
                // reset gradient to zero
                if upegged && sum < T::zero() {
                    for i in 0..self.m {
                        self.fjac[ij + i] = T::zero();
                    }
                }
            }
        }
    }

    ///
    /// On the first iteration and if user_scale is requested, scale according
    /// to the norms of the columns of the initial Jacobian,
    /// calculate the norm of the scaled x, and initialize the step bound delta.
    ///
    pub fn scale(&mut self) {
        if self.iter == 1 {
            if !self.cfg.do_user_scale {
                for j in 0..self.nfree {
                    self.diag[self.ifree[j]] = if self.wa2[j] == T::zero() {
                        T::one()
                    } else {
                        self.wa2[j]
                    };
                }
            }
            for j in 0..self.nfree {
                self.wa3[j] = self.diag[self.ifree[j]] * self.x[j];
            }
            self.xnorm = self.wa3.enorm();
            self.delta = self.cfg.step_factor * self.xnorm;
            if self.delta == T::zero() {
                self.delta = self.cfg.step_factor;
            }
        }
    }

    pub fn fill_xnew(&mut self) {
        for i in 0..self.nfree {
            self.xnew[self.ifree[i]] = self.x[i];
        }
    }

    ///
    /// Form (q transpose)*fvec and store the first n components in qtf.
    ///
    pub fn transpose(&mut self) {
        self.wa4.copy_from_slice(&self.fvec);
        let mut jj = 0;
        for j in 0..self.nfree {
            let temp = self.fjac[jj];
            if temp != T::zero() {
                let mut sum = T::zero();
                let mut ij = jj;
                for i in j..self.m {
                    sum += self.fjac[ij] * self.wa4[i];
                    ij += 1;
                }
                let temp = -sum / temp;
                ij = jj;
                for i in j..self.m {
                    self.wa4[i] += self.fjac[ij] * temp;
                    ij += 1;
                }
            }
            self.fjac[jj] = self.wa1[j];
            jj += self.m + 1;
            self.qtf[j] = self.wa4[j];
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
    pub fn gnorm(&self) -> T {
        let mut gnorm = T::zero();
        if self.fnorm != T::zero() {
            let mut jj = 0;
            for j in 0..self.nfree {
                let l = self.ipvt[j];
                if self.wa2[l] != T::zero() {
                    let mut sum = T::zero();
                    let mut ij = jj;
                    for i in 0..=j {
                        sum += self.fjac[ij] * (self.qtf[i] / self.fnorm);
                        ij += 1;
                    }
                    gnorm = gnorm.max((sum / self.wa2[l]).abs());
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
        let n_pegged = match self.f.model_params() {
            None => 0,
            Some(params) => {
                let mut n_pegged = 0;
                for (i, p) in params.iter().enumerate() {
                    if p.limit_lower.is_some_and(|x| x == self.xall[i])
                        || p.limit_upper.is_some_and(|x| x == self.xall[i])
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
            let l = j * self.m;
            for i in 0..self.nfree {
                covar[k + self.ifree[i]] = self.fjac[l + i]
            }
        }
        let mut xerror = vec![T::zero(); self.npar];
        for j in 0..self.nfree {
            let cc = self.fjac[j * self.m + j];
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
        let tolr = self.cfg.covtol * self.fjac[0].abs();
        let mut l: isize = -1;
        for k in 0..self.nfree {
            let k0 = k * self.m;
            let kk = k0 + k;
            if self.fjac[kk].abs() <= tolr {
                break;
            }
            self.fjac[kk] = T::one() / self.fjac[kk];
            for j in 0..k {
                let kj = k0 + j;
                let temp = self.fjac[kk] * self.fjac[kj];
                self.fjac[kj] = T::zero();
                let j0 = j * self.m;
                for i in 0..=j {
                    let dfjac = -temp * self.fjac[j0 + i];
                    self.fjac[k0 + i] += dfjac;
                }
            }
            l = k as isize;
        }
        // Form the full upper triangle of the inverse of (r transpose)*r
        // in the full upper triangle of r
        if l >= 0 {
            let l = l as usize;
            for k in 0..=l {
                let k0 = k * self.m;
                for j in 0..k {
                    let temp = self.fjac[k0 + j];
                    let j0 = j * self.m;
                    for i in 0..=j {
                        let dfjac = temp * self.fjac[k0 + i];
                        self.fjac[j0 + i] += dfjac;
                    }
                }
                let temp = self.fjac[k0 + k];
                for i in 0..=k {
                    self.fjac[k0 + i] *= temp;
                }
            }
        }
        // For the full lower triangle of the covariance matrix
        // in the strict lower triangle or and in wa
        for j in 0..self.nfree {
            let jj = self.ipvt[j];
            let sing = j as isize > l;
            let j0 = j * self.m;
            let jj0 = jj * self.m;
            for i in 0..=j {
                let ji = j0 + i;
                if sing {
                    self.fjac[ji] = T::zero();
                }
                let ii = self.ipvt[i];
                if ii > jj {
                    self.fjac[jj0 + ii] = self.fjac[ji];
                }
                if ii < jj {
                    self.fjac[ii * self.m + jj] = self.fjac[ji];
                }
            }
            self.wa2[jj] = self.fjac[j0 + j];
        }
        // Symmetrize the covariance matrix in r
        for j in 0..self.nfree {
            let j0 = j * self.m;
            for i in 0..j {
                self.fjac[j0 + i] = self.fjac[i * self.m + j];
            }
            self.fjac[j0 + j] = self.wa2[j];
        }
        self
    }

    pub fn rescale(&mut self) {
        if !self.cfg.do_user_scale {
            for j in 0..self.nfree {
                let i = self.ifree[j];
                self.diag[i] = self.diag[i].max(self.wa2[j]);
            }
        }
    }

    ///
    /// Given an m by nfree matrix a, an nfree by nfree nonsingular diagonal
    /// matrix d, an m-vector b, and a positive number delta,
    /// the problem is to determine a value for the parameter
    /// par such that if wa1 solves the system
    /// ```text
    ///     a*wa1 = b ,   sqrt(par)*d*wa1 = 0
    /// ```
    /// in the least squares sense, and dxnorm is the Euclidean
    /// norm of d*wa1, then either par is zero and
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
    pub fn lmpar(&mut self) {
        // Compute and store in wa1 the Gauss-Newton direction. If the
        // Jacobian is rank-deficient, obtain a least squares solution.
        let mut nsing = self.nfree;
        let mut jj = 0;
        for j in 0..self.nfree {
            self.wa3[j] = self.qtf[j];
            if self.fjac[jj] == T::zero() && nsing == self.nfree {
                nsing = j;
            }
            if nsing < self.nfree {
                self.wa3[j] = T::zero();
            }
            jj += self.m + 1;
        }
        if nsing >= 1 {
            for k in 0..nsing {
                let j = nsing - k - 1;
                let mut ij = self.m * j;
                self.wa3[j] /= self.fjac[j + ij];
                let temp = self.wa3[j];
                if j > 0 {
                    for i in 0..j {
                        self.wa3[i] -= self.fjac[ij] * temp;
                        ij += 1;
                    }
                }
            }
        }
        for j in 0..self.nfree {
            self.wa1[self.ipvt[j]] = self.wa3[j];
        }
        // Initialize the iteration counter.
        // Evaluate the function at the origin, and test
        // for acceptance of the Gauss-Newton direction.
        for j in 0..self.nfree {
            self.wa4[j] = self.diag[self.ifree[j]] * self.wa1[j];
        }
        let mut dxnorm = self.wa4[0..self.nfree].enorm();
        let mut fp = dxnorm - self.delta;
        if fp <= self.delta * T::P1 {
            self.par = T::zero();
            return;
        }
        // If the Jacobian is not rank deficient, the Newton
        // step provides a lower bound, parl, for the zero of
        // the function. Otherwise set this bound to zero.
        let mut parl = T::zero();
        if nsing >= self.nfree {
            self.newton_correction(dxnorm);
            let mut jj = 0;
            for j in 0..self.nfree {
                let mut sum = T::zero();
                if j > 0 {
                    let mut ij = jj;
                    for i in 0..j {
                        sum += self.fjac[ij] * self.wa3[i];
                        ij += 1;
                    }
                }
                self.wa3[j] = (self.wa3[j] - sum) / self.fjac[j + self.m * j];
                jj += self.m;
            }
            let temp = self.wa3[0..self.nfree].enorm();
            parl = ((fp / self.delta) / temp) / temp;
        }
        // Calculate an upper bound, paru, for the zero of the function.
        let mut jj = 0;
        for j in 0..self.nfree {
            let mut sum = T::zero();
            let mut ij = jj;
            for i in 0..=j {
                sum += self.fjac[ij] * self.qtf[i];
                ij += 1;
            }
            let l = self.ipvt[j];
            self.wa3[j] = sum / self.diag[self.ifree[l]];
            jj += self.m;
        }
        let gnorm = self.wa3[0..self.nfree].enorm();
        let mut paru = gnorm / self.delta;
        if paru == T::zero() {
            paru = T::MIN_POSITIVE / self.delta.min(T::P1);
        }
        // If the input par lies outside of the interval (parl,paru),
        // set par to the closer endpoint.
        self.par = self.par.max(parl);
        self.par = self.par.max(paru);
        if self.par == T::zero() {
            self.par = gnorm / dxnorm;
        }
        let mut iter = 0;
        loop {
            iter += 1;
            if self.par == T::zero() {
                self.par = T::MIN_POSITIVE.max(paru * T::P0001);
            }
            let temp = self.par.sqrt();
            for j in 0..self.nfree {
                self.wa3[j] = temp * self.diag[self.ifree[j]];
            }
            self.qrsolv();
            for j in 0..self.nfree {
                self.wa4[j] = self.diag[self.ifree[j]] * self.wa1[j];
            }
            dxnorm = self.wa4[0..self.nfree].enorm();
            let temp = fp;
            fp = dxnorm - self.delta;
            // If the function is small enough, accept the current value
            // of par. Also test for the exceptional cases where parl
            // is zero or the number of iterations has reached 10.
            if fp.abs() <= self.delta * T::P1
                || (parl == T::zero() && fp <= temp && temp < T::zero())
                || iter >= 10
            {
                return;
            }
            self.newton_correction(dxnorm);
            jj = 0;
            for j in 0..self.nfree {
                self.wa3[j] /= self.wa2[j];
                let temp = self.wa3[j];
                let jp1 = j + 1;
                if jp1 < self.nfree {
                    let mut ij = jp1 + jj;
                    for i in jp1..self.nfree {
                        self.wa3[i] -= self.fjac[ij] * temp;
                        ij += 1;
                    }
                }
                jj += self.m;
            }
            let temp = self.wa3[0..self.nfree].enorm();
            let parc = ((fp / self.delta) / temp) / temp;
            // Depending on the sign of the function, update parl or paru.
            if fp > T::zero() {
                parl = parl.max(self.par);
            }
            if fp < T::zero() {
                paru = paru.min(self.par);
            }
            // Compute an improved estimate for par.
            self.par = parl.max(self.par + parc);
        }
    }

    ///
    /// Compute the Newton correction.
    ///
    pub fn newton_correction(&mut self, dxnorm: T) {
        for j in 0..self.nfree {
            let l = self.ipvt[j];
            self.wa3[j] = self.diag[self.ifree[l]] * (self.wa4[l] / dxnorm);
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
    pub fn qrsolv(&mut self) {
        // Copy r and (q transpose)*b to preserve input and initialize s.
        // in particular, save the diagonal elements of r in x.
        let mut kk = 0;
        for j in 0..self.nfree {
            let mut ij = kk;
            let mut ik = kk;
            for _ in j..self.nfree {
                self.fjac[ij] = self.fjac[ik];
                ij += 1;
                ik += self.m;
            }
            self.wa1[j] = self.fjac[kk];
            self.wa4[j] = self.qtf[j];
            kk += self.m + 1;
        }
        // Eliminate the diagonal matrix d using a givens rotation.
        for j in 0..self.nfree {
            // Prepare the row of d to be eliminated, locating the
            // diagonal element using p from the qr factorization.
            let l = self.ipvt[j];
            if self.wa3[l] != T::zero() {
                for k in j..self.nfree {
                    self.wa2[k] = T::zero();
                }
                self.wa2[j] = self.wa3[l];
                // The transformations to eliminate the row of d
                // modify only a single element of (q transpose)*b
                // beyond the first n, which is initially zero.
                let mut qtbpj = T::zero();
                for k in j..self.nfree {
                    // Determine a givens rotation which eliminates the
                    // appropriate element in the current row of d.
                    if self.wa2[k] == T::zero() {
                        continue;
                    }
                    let kk = k + self.m * k;
                    let (sinx, cosx) = if self.fjac[kk].abs() < self.wa2[k].abs() {
                        let cotan = self.fjac[kk] / self.wa2[k];
                        let sinx = T::HALF / (cotan * cotan * T::P25 + T::P25).sqrt();
                        let cosx = sinx * cotan;
                        (sinx, cosx)
                    } else {
                        let tanx = self.wa2[k] / self.fjac[kk];
                        let cosx = T::HALF / (tanx * tanx * T::P25 + T::P25).sqrt();
                        let sinx = cosx * tanx;
                        (sinx, cosx)
                    };
                    // Compute the modified diagonal element of r and
                    // the modified element of ((q transpose)*b,0).
                    self.fjac[kk] = cosx * self.fjac[kk] + sinx * self.wa2[k];
                    let temp = cosx * self.wa4[k] + sinx * qtbpj;
                    qtbpj = -sinx * self.wa4[k] + cosx * qtbpj;
                    self.wa4[k] = temp;
                    // Accumulate the transformation in the row of s.
                    let kp1 = k + 1;
                    if self.nfree > kp1 {
                        let mut ik = kk + 1;
                        for i in kp1..self.nfree {
                            let temp = cosx * self.fjac[ik] + sinx * self.wa2[i];
                            self.wa2[i] = -sinx * self.fjac[ik] + cosx * self.wa2[i];
                            self.fjac[ik] = temp;
                            ik += 1;
                        }
                    }
                }
            }
            // Store the diagonal element of s and restore
            // the corresponding diagonal element of r.
            let kk = j + self.m * j;
            self.wa2[j] = self.fjac[kk];
            self.fjac[kk] = self.wa1[j];
        }
        // Solve the triangular system for z. if the system is
        // singular, then obtain a least squares solution.
        let mut nsing = self.nfree;
        for j in 0..self.nfree {
            if self.wa2[j] == T::zero() && nsing == self.nfree {
                nsing = j;
            }
            if nsing < self.nfree {
                self.wa4[j] = T::zero();
            }
        }
        if nsing > 0 {
            for k in 0..nsing {
                let j = nsing - k - 1;
                let mut sum = T::zero();
                let jp1 = j + 1;
                if nsing > jp1 {
                    let mut ij = jp1 + self.m * j;
                    for i in jp1..nsing {
                        sum += self.fjac[ij] * self.wa4[i];
                        ij += 1;
                    }
                }
                self.wa4[j] = (self.wa4[j] - sum) / self.wa2[j];
            }
        }
        // Permute the components of z back to components of x.
        for j in 0..self.nfree {
            self.wa1[self.ipvt[j]] = self.wa4[j];
        }
    }

    pub fn iterate(&mut self, gnorm: T) -> Result<MPFitDone, MPFitError> {
        for j in 0..self.nfree {
            self.wa1[j] = -self.wa1[j];
        }
        let mut alpha = T::one();
        if !self.qanylim {
            // No parameter limits, so just move to new position WA2
            for j in 0..self.nfree {
                self.wa2[j] = self.x[j] + self.wa1[j];
            }
        } else {
            // Respect the limits. If a step were to go out of bounds, then
            // we should take a step in the same direction but shorter distance.
            // The step should take us right to the limit in that case.
            for j in 0..self.nfree {
                let lpegged = self.qllim[j] && self.x[j] <= self.llim[j];
                let upegged = self.qulim[j] && self.x[j] >= self.ulim[j];
                let dwa1 = self.wa1[j].abs() > T::EPSILON;
                if lpegged && self.wa1[j] < T::zero() {
                    self.wa1[j] = T::zero();
                }
                if upegged && self.wa1[j] > T::zero() {
                    self.wa1[j] = T::zero();
                }
                if dwa1 && self.qllim[j] && self.x[j] + self.wa1[j] < self.llim[j] {
                    alpha = alpha.min((self.llim[j] - self.x[j]) / self.wa1[j]);
                }
                if dwa1 && self.qulim[j] && self.x[j] + self.wa1[j] > self.ulim[j] {
                    alpha = alpha.min((self.ulim[j] - self.x[j]) / self.wa1[j]);
                }
            }
            // Scale the resulting vector, advance to the next position
            for j in 0..self.nfree {
                self.wa1[j] *= alpha;
                self.wa2[j] = self.x[j] + self.wa1[j];
                // Adjust the output values. If the step put us exactly
                // on a boundary, make sure it is exact.
                let sgnu = if self.ulim[j] >= T::zero() {
                    T::one()
                } else {
                    -T::one()
                };
                let sgnl = if self.llim[j] >= T::zero() {
                    T::one()
                } else {
                    -T::one()
                };
                let ulim1 = self.ulim[j] * (T::one() - sgnu * T::EPSILON)
                    - if self.ulim[j] == T::zero() {
                        T::EPSILON
                    } else {
                        T::zero()
                    };
                let llim1 = self.llim[j] * (T::one() + sgnl * T::EPSILON)
                    + if self.llim[j] == T::zero() {
                        T::EPSILON
                    } else {
                        T::zero()
                    };
                if self.qulim[j] && self.wa2[j] >= ulim1 {
                    self.wa2[j] = self.ulim[j];
                }
                if self.qllim[j] && self.wa2[j] <= llim1 {
                    self.wa2[j] = self.llim[j];
                }
            }
        }
        for j in 0..self.nfree {
            self.wa3[j] = self.diag[self.ifree[j]] * self.wa1[j];
        }
        let pnorm = self.wa3[0..self.nfree].enorm();
        // On the first iteration, adjust the initial step bound.
        if self.iter == 1 {
            self.delta = self.delta.min(pnorm);
        }
        // Evaluate the function at x + p and calculate its norm.
        for i in 0..self.nfree {
            self.xnew[self.ifree[i]] = self.wa2[i];
        }
        self.f.evaluate(&self.xnew, &mut self.wa4)?;
        self.nfev += 1;
        self.fnorm1 = self.wa4[0..self.m].enorm();
        // Compute the scaled actual reduction.
        let actred = if self.fnorm1 * T::P1 < self.fnorm {
            let temp = self.fnorm1 / self.fnorm;
            T::one() - temp * temp
        } else {
            -T::one()
        };
        // Compute the scaled predicted reduction and
        // the scaled directional derivative.
        let mut jj = 0;
        for j in 0..self.nfree {
            self.wa3[j] = T::zero();
            let l = self.ipvt[j];
            let temp = self.wa1[l];
            let mut ij = jj;
            for i in 0..=j {
                self.wa3[i] += self.fjac[ij] * temp;
                ij += 1;
            }
            jj += self.m;
        }
        // Remember, alpha is the fraction of the full LM step actually
        // taken.
        let temp1 = self.wa3[0..self.nfree].enorm() * alpha / self.fnorm;
        let temp2 = ((alpha * self.par).sqrt() * pnorm) / self.fnorm;
        let temp11 = temp1 * temp1;
        let temp22 = temp2 * temp2;
        let prered = temp11 + temp22 / T::HALF;
        let dirder = -(temp11 + temp22);
        // Compute the ratio of the actual to the predicted
        // reduction.
        let ratio = if prered != T::zero() {
            actred / prered
        } else {
            T::zero()
        };
        // Update the step bound.
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
            self.par /= temp;
        } else if self.par == T::zero() || ratio >= T::P75 {
            self.delta = pnorm / T::HALF;
            self.par *= T::HALF;
        }
        // Test for successful iteration.
        if ratio >= T::P0001 {
            // Successful iteration. update x, fvec, and their norms.
            for j in 0..self.nfree {
                self.x[j] = self.wa2[j];
                self.wa2[j] = self.diag[self.ifree[j]] * self.x[j];
            }
            for i in 0..self.m {
                self.fvec[i] = self.wa4[i];
            }
            self.xnorm = self.wa2[0..self.nfree].enorm();
            self.fnorm = self.fnorm1;
            self.iter += 1;
        }
        // Tests for convergence.
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
        if self.info != MPFitInfo::NotDone {
            return Ok(MPFitDone::Exit);
        }
        // Tests for termination and stringent tolerances.
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
        if self.info != MPFitInfo::NotDone {
            return Ok(MPFitDone::Exit);
        }
        if ratio < T::P0001 {
            Ok(MPFitDone::Inner)
        } else {
            Ok(MPFitDone::Outer)
        }
    }

    pub fn check_config(&self) -> Result<(), MPFitError> {
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
}
