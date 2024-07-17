//!
//! Copied and modified from rmpfit crate Copyright (c) Vadim Dyadkin
//! Rust implementation of [CMPFIT](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html)
//!
use crate::fit::enums::MPFitSuccess;

#[derive(Debug, PartialEq, Eq)]
pub struct MPFitStatus<T> {
    /// Success enum
    pub success: MPFitSuccess,
    /// Final chi-squared
    pub best_norm: T,
    /// Starting value of chi-squared
    pub orig_norm: T,
    /// Number of iterations
    pub n_iter: usize,
    /// Number of function evaluations
    pub n_fev: usize,
    /// Total number of parameters
    pub n_par: usize,
    /// Number of free parameters
    pub n_free: usize,
    /// Number of pegged parameters
    pub n_pegged: usize,
    /// Number of residuals (= num. of data points)
    pub n_func: usize,
    /// Final residuals nfunc-vector
    pub residuals: Vec<T>,
    /// Final parameter uncertainties (1-sigma) npar-vector
    pub xerror: Vec<T>,
    /// Final parameter covariance matrix npar x npar array
    pub covar: Vec<T>,
}
