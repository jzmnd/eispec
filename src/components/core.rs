use num::traits::{ConstOne, ConstZero};

use crate::constants::FloatConst;
use crate::newtypes::{Current, Impedance, Voltage};

pub trait Component<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T>;
    fn calc_current(&self, v: Voltage<T>, freq: T) -> Current<T> {
        v / self.impedance(freq)
    }
    fn calc_voltage(&self, i: Current<T>, freq: T) -> Voltage<T> {
        i * self.impedance(freq)
    }
}
