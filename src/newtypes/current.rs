use num::complex::Complex;
use num::{Float, Num};
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

use crate::newtypes::impedance::Impedance;
use crate::newtypes::power::Power;
use crate::newtypes::voltage::Voltage;

///
/// Current (A)
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Current<T>(pub(super) Complex<T>);

impl<T> Current<T> {
    pub fn new(re: T, im: T) -> Self {
        Current(Complex::new(re, im))
    }
}

impl<T> fmt::Display for Current<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "I({}, {}j)", self.0.re, self.0.im)
    }
}

impl<T: Float> Current<T> {
    pub fn from_polar(r: T, theta: T) -> Self {
        Current(Complex::from_polar(r, theta))
    }
}

// Basic numerical derives
impl<T: Clone + Num> Add for Current<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl<T: Clone + Num> Sub for Current<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<T: Clone + Num> Mul for Current<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl<T: Clone + Num> Div for Current<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Power = Current * Voltage
impl<T: Clone + Num> Mul<Voltage<T>> for Current<T> {
    type Output = Power<T>;
    fn mul(self, other: Voltage<T>) -> Power<T> {
        Power::<T>(self.0 * other.0)
    }
}

// Voltage = Current * Impedance
impl<T: Clone + Num> Mul<Impedance<T>> for Current<T> {
    type Output = Voltage<T>;
    fn mul(self, other: Impedance<T>) -> Voltage<T> {
        Voltage::<T>(self.0 * other.0)
    }
}