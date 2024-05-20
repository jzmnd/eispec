use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::{Frequency, Impedance};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Cpe<T> {
    pub q0: T,
    pub n: T,
}

impl<T> Cpe<T> {
    pub fn new(q0: T, n: T) -> Self {
        Self { q0, n }
    }
}

impl<T> Component<T> for Cpe<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        let j = Complex::<T>::I;
        let z = Complex::from(self.q0) * (j * freq.to_angular()).powf(self.n).finv();
        Impedance::new(z.re, z.im)
    }
}
