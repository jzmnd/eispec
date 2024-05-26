use crate::constants::FloatConst;
use crate::newtypes::{Current, Frequency, Impedance, Voltage};

pub trait Component<T>
where
    T: FloatConst,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T>;
    fn calc_current(&self, v: Voltage<T>, freq: Frequency<T>) -> Current<T> {
        v / self.impedance(freq)
    }
    fn calc_voltage(&self, i: Current<T>, freq: Frequency<T>) -> Voltage<T> {
        i * self.impedance(freq)
    }
}
