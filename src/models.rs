use crate::constants::FloatConst;
use crate::newtypes::impedance::Impedance;
use num::complex::Complex;
use num::traits::{ConstOne, ConstZero};

pub fn freq_to_angular<T: FloatConst>(freq: T) -> T {
    T::PI_2 * freq
}

pub fn cpemodel<T>(r0: T, rinf: T, t0: T, alpha: T, freq: T) -> Impedance<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    let omega = freq_to_angular(freq);
    let one = Complex::<T>::ONE;
    let j = Complex::<T>::I;
    let b = one + (j * omega * t0).powf(T::one() - alpha);
    let zc = Complex::from(r0 - rinf) / b;
    let z = Complex::from(rinf) + zc;
    Impedance::new(z.re, z.im)
}

pub fn warburg_inf<T>(r0: T, rinf: T, t0: T, freq: T) -> Impedance<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    let omega = freq_to_angular(freq);
    let one = Complex::<T>::ONE;
    let j = Complex::<T>::I;
    let b = one + (j * omega * t0).sqrt();
    let zc = Complex::from(r0 - rinf) / b;
    let z = Complex::from(rinf) + zc;
    Impedance::new(z.re, z.im)
}

pub fn cpemodel_beta<T>(r0: T, rinf: T, t0: T, alpha: T, beta: T, freq: T) -> Impedance<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    let omega = freq_to_angular(freq);
    let one = Complex::<T>::ONE;
    let j = Complex::<T>::I;
    let b = (one + (j * omega * t0).powf(T::one() - alpha)).powf(beta);
    let zc = Complex::from(r0 - rinf) / b;
    let z = Complex::from(rinf) + zc;
    Impedance::new(z.re, z.im)
}

pub fn hnmodel<T>(r0: T, rinf: T, t0: T, alpha: T, beta: T, freq: T) -> Impedance<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    let omega = freq_to_angular(freq);
    let one = Complex::<T>::ONE;
    let j = Complex::<T>::I;
    let b = (one + (j * omega * t0).powf(T::one() - alpha)).powf(beta);
    let zc = Complex::from(r0 - rinf) / (Complex::from(r0 / rinf) * (b - one) + one);
    let z = Complex::from(rinf) + zc;
    Impedance::new(z.re, z.im)
}

pub fn cpemodel_interlayer<T>(r0: T, rinf: T, cinter: T, t0: T, alpha: T, freq: T) -> Impedance<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    let omega = freq_to_angular(freq);
    let one = Complex::<T>::ONE;
    let j = Complex::<T>::I;
    let b = one + (j * omega * t0).powf(T::one() - alpha);
    let zc = Complex::from(r0 - rinf) / b;
    let zinter = one / (j * omega * cinter);
    let z = Complex::from(rinf) + zc + zinter;
    Impedance::new(z.re, z.im)
}

pub fn oxidemodel<T>(rs: T, q1: T, n1: T, r1: T, q2: T, n2: T, r2: T, freq: T) -> Impedance<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    let omega = freq_to_angular(freq);
    let one = Complex::<T>::ONE;
    let j = Complex::<T>::I;
    let cperp1 = Complex::from(r1) / (Complex::from(r1) * (j * omega * q1).powf(n1) + one);
    let cperp2 = Complex::from(r2) / (Complex::from(r2) * (j * omega * q2).powf(n2) + one);
    let z = Complex::from(rs) + cperp1 + cperp2;
    Impedance::new(z.re, z.im)
}
