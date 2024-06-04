//!
//! Floating point constants.
//!
//! Some of these are specific to the CMPFIT algorithm so are defined
//! here rather than just using `num::traits::FloatConst`.
//!
use num::traits::{ConstOne, ConstZero};
use num::Float;

pub trait FloatConst: Float + ConstOne + ConstZero {
    const PI: Self;
    const PI_2: Self;
    const SQRT_2: Self;
    const FRAC_2_PI: Self;
    const FRAC_PI_2: Self;
    const EPSILON: Self;
    const MP_RDWARF: Self;
    const MP_RGIANT: Self;
    const MIN_POSITIVE: Self;
    const HALF: Self;
    const P1: Self;
    const P25: Self;
    const P75: Self;
    const P05: Self;
    const P0001: Self;
}

impl FloatConst for f32 {
    const PI: f32 = std::f32::consts::PI;
    const PI_2: f32 = std::f32::consts::TAU;
    const SQRT_2: f32 = std::f32::consts::SQRT_2;
    const FRAC_2_PI: f32 = std::f32::consts::FRAC_2_PI;
    const FRAC_PI_2: f32 = std::f32::consts::FRAC_PI_2;
    const EPSILON: f32 = std::f32::EPSILON;
    const MP_RDWARF: f32 = 1.3278711e-18;
    const MP_RGIANT: f32 = 1.8446743e18;
    const MIN_POSITIVE: f32 = std::f32::MIN_POSITIVE;
    const HALF: f32 = 0.5;
    const P1: f32 = 0.1;
    const P25: f32 = 0.25;
    const P75: f32 = 0.75;
    const P05: f32 = 0.05;
    const P0001: f32 = 1e-4;
}

impl FloatConst for f64 {
    const PI: f64 = std::f64::consts::PI;
    const PI_2: f64 = std::f64::consts::TAU;
    const SQRT_2: f64 = std::f64::consts::SQRT_2;
    const FRAC_2_PI: f64 = std::f64::consts::FRAC_2_PI;
    const FRAC_PI_2: f64 = std::f64::consts::FRAC_PI_2;
    const EPSILON: f64 = std::f64::EPSILON;
    const MP_RDWARF: f64 = 1.8269129119256895e-153;
    const MP_RGIANT: f64 = 1.3407807929942596e153;
    const MIN_POSITIVE: f64 = std::f64::MIN_POSITIVE;
    const HALF: f64 = 0.5;
    const P1: f64 = 0.1;
    const P25: f64 = 0.25;
    const P75: f64 = 0.75;
    const P05: f64 = 0.05;
    const P0001: f64 = 1e-4;
}
