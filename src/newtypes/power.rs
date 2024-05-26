use num::complex::Complex;
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, Div, Mul, Sub};

use crate::constants::FloatConst;
use crate::newtypes::{Current, Voltage};

///
/// Power (W)
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Power<T>(pub(super) Complex<T>);

impl<T> Power<T> {
    pub fn new(re: T, im: T) -> Self {
        Self(Complex::new(re, im))
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

impl<T: FloatConst> Power<T> {
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
impl<T: Clone + FloatConst> Add for Power<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl<T: Clone + FloatConst> Sub for Power<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<T: Clone + FloatConst> Mul for Power<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl<T: Clone + FloatConst> Div for Power<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Sum derives
impl<T: Clone + FloatConst> Sum for Power<T> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let pzero = Power::new(T::zero(), T::zero());
        iter.fold(pzero, |acc, c| acc + c)
    }
}

// Current = Power / Voltage
impl<T: Clone + FloatConst> Div<Voltage<T>> for Power<T> {
    type Output = Current<T>;
    fn div(self, other: Voltage<T>) -> Current<T> {
        Current::<T>(self.0 / other.0)
    }
}

// Voltage = Power / Current
impl<T: Clone + FloatConst> Div<Current<T>> for Power<T> {
    type Output = Voltage<T>;
    fn div(self, other: Current<T>) -> Voltage<T> {
        Voltage::<T>(self.0 / other.0)
    }
}
