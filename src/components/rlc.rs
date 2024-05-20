use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::{Frequency, Impedance};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Resistor<T> {
    pub r0: T,
}

impl<T> Resistor<T> {
    pub fn new(r0: T) -> Self {
        Self { r0 }
    }
}

impl<T> Component<T> for Resistor<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, _freq: Frequency<T>) -> Impedance<T> {
        Impedance::new(self.r0, T::zero())
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Capacitor<T> {
    pub c0: T,
}

impl<T> Capacitor<T> {
    pub fn new(c0: T) -> Self {
        Self { c0 }
    }
}

impl<T> Component<T> for Capacitor<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        Impedance::new(T::zero(), -(freq.to_angular() * self.c0).recip())
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct Inductor<T> {
    pub l0: T,
}

impl<T> Inductor<T> {
    pub fn new(l0: T) -> Self {
        Self { l0 }
    }
}

impl<T> Component<T> for Inductor<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        Impedance::new(T::zero(), freq.to_angular() * self.l0)
    }
}
