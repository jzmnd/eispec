//!
//! Copied and modified from rmpfit crate Copyright (c) Vadim Dyadkin
//! Rust implementation of [CMPFIT](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html)
//!
use crate::constants::FloatConst;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct MPFitConfig<T> {
    /// Relative chi-square convergence criterion
    pub ftol: T,
    /// Relative parameter convergence criterion
    pub xtol: T,
    /// Orthogonality convergence criterion
    pub gtol: T,
    /// Finite derivative step size
    pub epsfcn: T,
    /// Initial step bound
    pub step_factor: T,
    /// Range tolerance for covariance calculation
    pub covtol: T,
    /// Maximum number of iterations
    /// If max_iter == 0, then basic error checking is done, and
    /// parameter errors/covariances are estimated based on input
    /// parameter values, but no fitting iterations are done.
    pub max_iter: usize,
    /// Maximum number of function evaluations
    /// If max_fev == 0 no limit is applied
    pub max_fev: usize,
    /// Scale variables by user values?
    /// true = yes, user scale values in diag;
    /// false = no, variables scaled internally
    pub do_user_scale: bool,
    /// Disable check for infinite quantities from user?
    /// true = perform check;
    /// false = do not perform check
    pub no_finite_check: bool,
}

impl<T> Default for MPFitConfig<T>
where
    T: FloatConst,
{
    fn default() -> Self {
        Self {
            ftol: T::from(1e-10).unwrap(),
            xtol: T::from(1e-10).unwrap(),
            gtol: T::from(1e-10).unwrap(),
            epsfcn: T::EPSILON,
            step_factor: T::from(100.0).unwrap(),
            covtol: T::from(1e-14).unwrap(),
            max_iter: 200,
            max_fev: 0,
            do_user_scale: false,
            no_finite_check: false,
        }
    }
}
