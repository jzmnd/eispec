use num::complex::Complex;
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, Div, Mul, Sub};

use crate::constants::FloatConst;
use crate::newtypes::{Current, Voltage};

///
/// Impedance (Ohm)
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Impedance<T>(pub(super) Complex<T>);

impl<T> Impedance<T> {
    pub fn new(re: T, im: T) -> Self {
        Self(Complex::new(re, im))
    }
}

impl<T> fmt::Display for Impedance<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Z({}, {}j)", self.0.re, self.0.im)
    }
}

impl<T: FloatConst> Impedance<T> {
    pub fn re(&self) -> T {
        self.0.re
    }

    pub fn im(&self) -> T {
        self.0.im
    }

    pub fn from_polar(r: T, theta: T) -> Self {
        Self(Complex::from_polar(r, theta))
    }

    pub fn to_polar(self) -> (T, T) {
        self.0.to_polar()
    }

    pub fn from_cpg(cp: T, g: T, freq: T) -> Self {
        let omega = T::PI_2 * freq;
        let r = g.recip();
        let re = r / (T::one() + (omega * r * cp).powi(2));
        let im = -omega * cp * (r.powi(2)) / (T::one() + (omega * r * cp).powi(2));
        Self::new(re, im)
    }

    pub fn from_csg(cs: T, g: T, freq: T) -> Self {
        let omega = T::PI_2 * freq;
        let r = g.recip();
        let im = -(omega * cs).recip();
        Self::new(r, im)
    }

    pub fn finv(self) -> Self {
        Self(self.0.finv())
    }
}

// Basic numerical derives
impl<T: Clone + FloatConst> Add for Impedance<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl<T: Clone + FloatConst> Sub for Impedance<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<T: Clone + FloatConst> Mul for Impedance<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl<T: Clone + FloatConst> Div for Impedance<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Sum derives
impl<T: Clone + FloatConst> Sum for Impedance<T> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let zzero = Impedance::new(T::zero(), T::zero());
        iter.fold(zzero, |acc, c| acc + c)
    }
}

// Voltage = Impedance * Current
impl<T: Clone + FloatConst> Mul<Current<T>> for Impedance<T> {
    type Output = Voltage<T>;
    fn mul(self, other: Current<T>) -> Voltage<T> {
        Voltage::<T>(self.0 * other.0)
    }
}
