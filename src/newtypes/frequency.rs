use num::{Float, Num};
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, Div, Mul, Sub};

use crate::constants::FloatConst;

///
/// Frequency (Hz)
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Frequency<T>(pub(super) T);

impl<T> Frequency<T> {
    pub fn new(f: T) -> Self {
        Self(f)
    }
}

impl<T> fmt::Display for Frequency<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Freq({})", self.0)
    }
}

impl<T: Float + FloatConst> Frequency<T> {
    pub fn from_angular(omega: T) -> Self {
        Self(omega / T::PI_2)
    }

    pub fn to_angular(self) -> T {
        self.0 * T::PI_2
    }

    pub fn value(self) -> T {
        self.0
    }
}

// Basic numerical derives
impl<T: Clone + Num> Add for Frequency<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl<T: Clone + Num> Sub for Frequency<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<T: Clone + Num> Mul for Frequency<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl<T: Clone + Num> Div for Frequency<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Sum derives
impl<T: Clone + Num> Sum for Frequency<T> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Frequency(T::zero()), |acc, c| acc + c)
    }
}
