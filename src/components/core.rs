//!
//! Core features of an electrical component or circuit element.
//!
//! The Component trait should be implemented for an circuit element such
//! as a resistor or a constant phase element. The only function that needs
//! to be defined is the frequency dependent complex impedance.
//!
use crate::constants::FloatConst;
use crate::newtypes::{Current, Frequency, Impedance, Voltage};

pub trait Component<T>
where
    T: FloatConst,
{
    /// Returns the complex impedance of the component at a given frequency
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T>;
    /// Returns the complex current through the component at a given frequency and voltage
    fn calc_current(&self, v: Voltage<T>, freq: Frequency<T>) -> Current<T> {
        v / self.impedance(freq)
    }
    /// Returns the complex voltage across the component at a given frequency and current
    fn calc_voltage(&self, i: Current<T>, freq: Frequency<T>) -> Voltage<T> {
        i * self.impedance(freq)
    }
}
