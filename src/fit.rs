//!
//! Complex impedance fitting using the Levenberg-Marquardt algorithm.
//!
pub mod config;
pub mod enorm;
pub mod enums;
pub mod jacobian;
mod model;
mod model_parameter;
pub mod mpfit;
pub mod status;

pub use crate::fit::model::ImpedanceModel;
pub use crate::fit::model_parameter::{ModelParameter, ParameterBounds};
