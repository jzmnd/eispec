use num::complex::Complex;
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, Div, Mul, Sub};

use crate::constants::FloatConst;
use crate::newtypes::{Impedance, Power, Voltage};

///
/// Current (A)
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Current<T>(pub(super) Complex<T>);

impl<T> Current<T> {
    pub fn new(re: T, im: T) -> Self {
        Self(Complex::new(re, im))
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

impl<T: FloatConst> Current<T> {
    pub fn re(&self) -> T {
        self.0.re
    }

    pub fn im(&self) -> T {
        self.0.im
    }

    pub fn from_polar(r: T, theta: T) -> Self {
        Self(Complex::from_polar(r, theta))
    }

    pub fn finv(self) -> Self {
        Self(self.0.finv())
    }
}

// Basic numerical derives
impl<T: Clone + FloatConst> Add for Current<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl<T: Clone + FloatConst> Sub for Current<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<T: Clone + FloatConst> Mul for Current<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl<T: Clone + FloatConst> Div for Current<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Sum derives
impl<T: Clone + FloatConst> Sum for Current<T> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let izero = Current::new(T::zero(), T::zero());
        iter.fold(izero, |acc, c| acc + c)
    }
}

// Power = Current * Voltage
impl<T: Clone + FloatConst> Mul<Voltage<T>> for Current<T> {
    type Output = Power<T>;
    fn mul(self, other: Voltage<T>) -> Power<T> {
        Power::<T>(self.0 * other.0)
    }
}

// Voltage = Current * Impedance
impl<T: Clone + FloatConst> Mul<Impedance<T>> for Current<T> {
    type Output = Voltage<T>;
    fn mul(self, other: Impedance<T>) -> Voltage<T> {
        Voltage::<T>(self.0 * other.0)
    }
}
