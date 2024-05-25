//!
//! Warburg Elements.
//!
use num::complex::Complex;

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::{Frequency, Impedance};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Warburg<T> {
    pub aw: T,
}

impl<T> Warburg<T> {
    pub fn new(aw: T) -> Self {
        Self { aw }
    }
}

impl<T> Component<T> for Warburg<T>
where
    T: FloatConst,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        let j = Complex::<T>::I;
        let sqrtjo = (j * freq.to_angular()).sqrt();
        let z = Complex::from(T::SQRT_2 * self.aw) / sqrtjo;
        Impedance::new(z.re, z.im)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct WarburgShort<T> {
    pub aw: T,
    pub d: T,
    pub diffusion_coeff: T,
}

impl<T> WarburgShort<T> {
    pub fn new(aw: T, d: T, diffusion_coeff: T) -> Self {
        Self {
            aw,
            d,
            diffusion_coeff,
        }
    }
}

impl<T> Component<T> for WarburgShort<T>
where
    T: FloatConst,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        let j = Complex::<T>::I;
        let sqrtjo = (j * freq.to_angular()).sqrt();
        let b = self.d / self.diffusion_coeff.sqrt();
        let z = Complex::from(self.aw) * (sqrtjo * b).tanh() / sqrtjo;
        Impedance::new(z.re, z.im)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct WarburgOpen<T> {
    pub aw: T,
    pub d: T,
    pub diffusion_coeff: T,
}

impl<T> WarburgOpen<T> {
    pub fn new(aw: T, d: T, diffusion_coeff: T) -> Self {
        Self {
            aw,
            d,
            diffusion_coeff,
        }
    }
}

impl<T> Component<T> for WarburgOpen<T>
where
    T: FloatConst,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        let j = Complex::<T>::I;
        let sqrtjo = (j * freq.to_angular()).sqrt();
        let b = self.d / self.diffusion_coeff.sqrt();
        let z = Complex::from(self.aw) / (sqrtjo * b).tanh() / sqrtjo;
        Impedance::new(z.re, z.im)
    }
}
