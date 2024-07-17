pub mod config;
pub mod enorm;
pub mod enums;
mod fit;
pub mod mpfit;
pub mod status;

pub use crate::fit::fit::{ImpedanceDataFitter, ModelParameter};
