use crate::constants::FloatConst;
use crate::newtypes::impedance::Impedance;
use crate::utils::freq_to_angular;
use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

pub trait Component<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T>;
}

pub struct Resistor<T> {
    pub r0: T,
}

impl<T> Component<T> for Resistor<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, _freq: T) -> Impedance<T> {
        Impedance::new(self.r0, T::zero())
    }
}

pub struct Capacitor<T> {
    pub c0: T,
}

impl<T> Component<T> for Capacitor<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        Impedance::new(T::zero(), (omega * self.c0).recip())
    }
}

pub struct Inductor<T> {
    pub l0: T,
}

impl<T> Component<T> for Inductor<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        let omega = freq_to_angular(freq);
        Impedance::new(T::zero(), omega * self.l0)
    }
}

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
        let one = Complex::<T>::ONE;
        let z = one / Complex::from(self.q0) * (j * omega).powf(self.n);
        Impedance::new(z.re, z.im)
    }
}

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
        let sqrttwo = T::from(2).unwrap().sqrt();
        let z = Complex::from(sqrttwo * self.aw) / (j * omega).sqrt();
        Impedance::new(z.re, z.im)
    }
}

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
        let b = self.d / self.diffusion_coeff.sqrt();
        let z = Complex::from(self.aw) * ((j * omega).sqrt() * b).tanh() / (j * omega).sqrt();
        Impedance::new(z.re, z.im)
    }
}

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
        let b = self.d / self.diffusion_coeff.sqrt();
        let z = Complex::from(self.aw) / ((j * omega).sqrt() * b).tanh() / (j * omega).sqrt();
        Impedance::new(z.re, z.im)
    }
}

pub struct Gerischer<T> {
    pub rg: T,
    pub tau: T,
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
