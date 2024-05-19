use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::Impedance;
use crate::utils::freq_to_angular;

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Warburg<T> {
    pub aw: T,
}

impl<T> Component<T> for Warburg<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let sqrtjo = (j * omega).sqrt();
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

impl<T> Component<T> for WarburgShort<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let sqrtjo = (j * omega).sqrt();
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

impl<T> Component<T> for WarburgOpen<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let sqrtjo = (j * omega).sqrt();
        let b = self.d / self.diffusion_coeff.sqrt();
        let z = Complex::from(self.aw) / (sqrtjo * b).tanh() / sqrtjo;
        Impedance::new(z.re, z.im)
    }
}
