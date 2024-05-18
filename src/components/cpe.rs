use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::Impedance;
use crate::utils::freq_to_angular;

pub struct Cpe<T> {
    pub q0: T,
    pub n: T,
}

impl<T> Component<T> for Cpe<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let z = Complex::from(self.q0) * (j * omega).powf(self.n).finv();
        Impedance::new(z.re, z.im)
    }
}
