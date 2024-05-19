use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::Impedance;
use crate::utils::freq_to_angular;

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Gerischer<T> {
    pub rg: T,
    pub tau: T,
}

impl<T> Gerischer<T> {
    pub fn new(rg: T, tau: T) -> Self {
        Self { rg, tau }
    }
}

impl<T> Component<T> for Gerischer<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let one = Complex::<T>::ONE;
        let z = Complex::from(self.rg) / (one + (j * omega * self.tau).sqrt());
        Impedance::new(z.re, z.im)
    }
}
