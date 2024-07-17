//!
//! Copied and modified from rmpfit crate Copyright (c) Vadim Dyadkin
//! Rust implementation of [CMPFIT](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html)
//!
use std::fmt;
use thiserror::Error;

#[derive(Debug, PartialEq, Eq)]
pub enum MPFitInfo {
    /// Not finished iterations
    NotDone,
    /// Convergence in chi-square value
    ConvergenceChi,
    /// Convergence in parameter value
    ConvergencePar,
    /// Convergence in both chi-square and parameter
    ConvergenceBoth,
    /// Convergence in orthogonality
    ConvergenceDir,
    /// Maximum number of iterations reached
    MaxIterReached,
    /// ftol is too small; no further improvement
    FtolNoImprovement,
    /// xtol is too small; no further improvement
    XtolNoImprovement,
    /// gtol is too small; no further improvement
    GtolNoImprovement,
}

impl fmt::Display for MPFitInfo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::NotDone => "Unknown error",
                Self::ConvergenceChi => "Convergence in chi-square value",
                Self::ConvergencePar => "Convergence in parameter value",
                Self::ConvergenceBoth => "Convergence in chi-square and parameter values",
                Self::ConvergenceDir => "Convergence in orthogonality",
                Self::MaxIterReached => "Maximum number of iterations reached",
                Self::FtolNoImprovement => "ftol is too small; no further improvement",
                Self::XtolNoImprovement => "xtol is too small; no further improvement",
                Self::GtolNoImprovement => "gtol is too small; no further improvement",
            }
        )
    }
}

#[derive(Debug, Error)]
pub enum MPFitError {
    #[error("General input parameter error")]
    Input,
    #[error("User function produced non-finite values")]
    Nan,
    #[error("No user data points were supplied")]
    Empty,
    #[error("No free parameters")]
    NoFree,
    #[error("Initial values inconsistent with constraints")]
    InitBounds,
    #[error("Initial constraints inconsistent")]
    Bounds,
    #[error("Not enough degrees of freedom")]
    DoF,
    #[error("Error during user evaluation")]
    Eval,
}

pub enum MPFitDone {
    /// Fit is complete, exit
    Exit,
    /// Fit is running the inner loop
    Inner,
    /// Fit is running the outer loop
    Outer,
}
