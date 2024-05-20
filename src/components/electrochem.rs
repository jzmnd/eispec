use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::{Frequency, Impedance};

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
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        let j = Complex::<T>::I;
        let one = Complex::<T>::ONE;
        let z = Complex::from(self.rg) / (one + (j * freq.to_angular() * self.tau).sqrt());
        Impedance::new(z.re, z.im)
    }
}
