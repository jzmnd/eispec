use num::Float;

pub trait FloatConst: Float {
    const PI: Self;
    const PI_2: Self;
    const SQRT_2: Self;
    const FRAC_2_PI: Self;
    const FRAC_PI_2: Self;
}

impl FloatConst for f32 {
    const PI: f32 = std::f32::consts::PI;
    const PI_2: f32 = std::f32::consts::TAU;
    const SQRT_2: f32 = std::f32::consts::SQRT_2;
    const FRAC_2_PI: f32 = std::f32::consts::FRAC_2_PI;
    const FRAC_PI_2: f32 = std::f32::consts::FRAC_PI_2;
}

impl FloatConst for f64 {
    const PI: f64 = std::f64::consts::PI;
    const PI_2: f64 = std::f64::consts::TAU;
    const SQRT_2: f64 = std::f64::consts::SQRT_2;
    const FRAC_2_PI: f64 = std::f64::consts::FRAC_2_PI;
    const FRAC_PI_2: f64 = std::f64::consts::FRAC_PI_2;
}
