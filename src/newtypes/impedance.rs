use num::complex::Complex;
use num::{Float, Num};
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

use crate::constants::FloatConst;
use crate::newtypes::current::Current;
use crate::newtypes::voltage::Voltage;

///
/// Impedance (Ohm)
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Impedance<T>(pub(super) Complex<T>);

impl<T> Impedance<T> {
    pub fn new(re: T, im: T) -> Self {
        Impedance(Complex::new(re, im))
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

impl<T: Float + FloatConst> Impedance<T> {
    pub fn from_polar(r: T, theta: T) -> Self {
        Impedance(Complex::from_polar(r, theta))
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
}

// Basic numerical derives
impl<T: Clone + Num> Add for Impedance<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl<T: Clone + Num> Sub for Impedance<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<T: Clone + Num> Mul for Impedance<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl<T: Clone + Num> Div for Impedance<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Voltage = Impedance * Current
impl<T: Clone + Num> Mul<Current<T>> for Impedance<T> {
    type Output = Voltage<T>;
    fn mul(self, other: Current<T>) -> Voltage<T> {
        Voltage::<T>(self.0 * other.0)
    }
}