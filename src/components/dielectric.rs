use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::Impedance;
use crate::utils::freq_to_angular;

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct HavriliakNegami<T> {
    pub rinf: T,
    pub r0: T,
    pub tau: T,
    pub alpha: T,
    pub beta: T,
}

impl<T> HavriliakNegami<T> {
    pub fn new(rinf: T, r0: T, tau: T, alpha: T, beta: T) -> Self {
        Self {
            rinf,
            r0,
            tau,
            alpha,
            beta,
        }
    }
}

impl<T> Component<T> for HavriliakNegami<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let jot = j * omega * self.tau;
        let one = Complex::<T>::ONE;
        let zc = Complex::from(self.r0) / (one + jot.powf(self.alpha)).powf(self.beta);
        let z = zc + self.rinf;
        Impedance::new(z.re, z.im)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Debye<T> {
    pub rinf: T,
    pub r0: T,
    pub tau: T,
}

impl<T> Debye<T> {
    pub fn new(rinf: T, r0: T, tau: T) -> Self {
        Self { rinf, r0, tau }
    }
}

impl<T> Component<T> for Debye<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let jot = j * omega * self.tau;
        let one = Complex::<T>::ONE;
        let zc = Complex::from(self.r0) / (one + jot);
        let z = zc + self.rinf;
        Impedance::new(z.re, z.im)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct ColeCole<T> {
    pub rinf: T,
    pub r0: T,
    pub tau: T,
    pub alpha: T,
}

impl<T> ColeCole<T> {
    pub fn new(rinf: T, r0: T, tau: T, alpha: T) -> Self {
        Self {
            rinf,
            r0,
            tau,
            alpha,
        }
    }
}

impl<T> Component<T> for ColeCole<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let jot = j * omega * self.tau;
        let one = Complex::<T>::ONE;
        let zc = Complex::from(self.r0) / (one + jot.powf(self.alpha));
        let z = zc + self.rinf;
        Impedance::new(z.re, z.im)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct ColeDavidson<T> {
    pub rinf: T,
    pub r0: T,
    pub tau: T,
    pub beta: T,
}

impl<T> ColeDavidson<T> {
    pub fn new(rinf: T, r0: T, tau: T, beta: T) -> Self {
        Self {
            rinf,
            r0,
            tau,
            beta,
        }
    }
}

impl<T> Component<T> for ColeDavidson<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        let j = Complex::<T>::I;
        let jot = j * omega * self.tau;
        let one = Complex::<T>::ONE;
        let zc = Complex::from(self.r0) / (one + jot).powf(self.beta);
        let z = zc + self.rinf;
        Impedance::new(z.re, z.im)
    }
}
