use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::Impedance;
use crate::utils::freq_to_angular;

pub struct HavriliakNegami<T> {
    pub rinf: T,
    pub r0: T,
    pub tau: T,
    pub alpha: T,
    pub beta: T,
}

impl<T> Component<T> for HavriliakNegami<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let one = Complex::<T>::ONE;
        let zc = Complex::from(self.r0)
            / (one + (j * omega * self.tau).powf(self.alpha)).powf(self.beta);
        let z = zc + self.rinf;
        Impedance::new(z.re, z.im)
    }
}

pub struct Debye<T> {
    pub rinf: T,
    pub r0: T,
    pub tau: T,
}

impl<T> Component<T> for Debye<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let one = Complex::<T>::ONE;
        let zc = Complex::from(self.r0) / (one + (j * omega * self.tau));
        let z = zc + self.rinf;
        Impedance::new(z.re, z.im)
    }
}

pub struct ColeCole<T> {
    pub rinf: T,
    pub r0: T,
    pub tau: T,
    pub alpha: T,
}

impl<T> Component<T> for ColeCole<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let one = Complex::<T>::ONE;
        let zc = Complex::from(self.r0) / (one + (j * omega * self.tau).powf(self.alpha));
        let z = zc + self.rinf;
        Impedance::new(z.re, z.im)
    }
}

pub struct ColeDavidson<T> {
    pub rinf: T,
    pub r0: T,
    pub tau: T,
    pub beta: T,
}

impl<T> Component<T> for ColeDavidson<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let one = Complex::<T>::ONE;
        let zc = Complex::from(self.r0) / (one + (j * omega * self.tau)).powf(self.beta);
        let z = zc + self.rinf;
        Impedance::new(z.re, z.im)
    }
}
