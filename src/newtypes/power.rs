use num::complex::Complex;
use num::{Float, Num};
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

use crate::newtypes::current::Current;
use crate::newtypes::voltage::Voltage;

///
/// Power (W)
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Power<T>(pub(super) Complex<T>);

impl<T> Power<T> {
    pub fn new(re: T, im: T) -> Self {
        Power(Complex::new(re, im))
    }
}

impl<T> fmt::Display for Power<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "P({}, {}j)", self.0.re, self.0.im)
    }
}

impl<T: Float> Power<T> {
    pub fn from_polar(r: T, theta: T) -> Self {
        Power(Complex::from_polar(r, theta))
    }
}

// Basic numerical derives
impl<T: Clone + Num> Add for Power<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl<T: Clone + Num> Sub for Power<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<T: Clone + Num> Mul for Power<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl<T: Clone + Num> Div for Power<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Current = Power / Voltage
impl<T: Clone + Num> Div<Voltage<T>> for Power<T> {
    type Output = Current<T>;
    fn div(self, other: Voltage<T>) -> Current<T> {
        Current::<T>(self.0 / other.0)
    }
}

// Voltage = Power / Current
impl<T: Clone + Num> Div<Current<T>> for Power<T> {
    type Output = Voltage<T>;
    fn div(self, other: Current<T>) -> Voltage<T> {
        Voltage::<T>(self.0 / other.0)
    }
}