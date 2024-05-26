use num::complex::Complex;
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, Div, Mul, Sub};

use crate::constants::FloatConst;
use crate::newtypes::{Current, Impedance, Power};

///
/// Voltage (V)
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Voltage<T>(pub(super) Complex<T>);

impl<T> Voltage<T> {
    pub fn new(re: T, im: T) -> Self {
        Self(Complex::new(re, im))
    }
}

impl<T> fmt::Display for Voltage<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "V({}, {}j)", self.0.re, self.0.im)
    }
}

impl<T: FloatConst> Voltage<T> {
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
impl<T: Clone + FloatConst> Add for Voltage<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl<T: Clone + FloatConst> Sub for Voltage<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<T: Clone + FloatConst> Mul for Voltage<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl<T: Clone + FloatConst> Div for Voltage<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Sum derives
impl<T: Clone + FloatConst> Sum for Voltage<T> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let vzero = Voltage::new(T::zero(), T::zero());
        iter.fold(vzero, |acc, c| acc + c)
    }
}

// Power = Voltage * Current
impl<T: Clone + FloatConst> Mul<Current<T>> for Voltage<T> {
    type Output = Power<T>;
    fn mul(self, other: Current<T>) -> Power<T> {
        Power::<T>(self.0 * other.0)
    }
}

// Current = Voltage / Impedance
impl<T: Clone + FloatConst> Div<Impedance<T>> for Voltage<T> {
    type Output = Current<T>;
    fn div(self, other: Impedance<T>) -> Current<T> {
        Current::<T>(self.0 / other.0)
    }
}

// Impedance = Voltage / Current
impl<T: Clone + FloatConst> Div<Current<T>> for Voltage<T> {
    type Output = Impedance<T>;
    fn div(self, other: Current<T>) -> Impedance<T> {
        Impedance::<T>(self.0 / other.0)
    }
}
