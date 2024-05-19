use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::Impedance;

#[derive(Default)]
pub struct ParallelCircuit<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    pub components: Vec<Box<dyn Component<T>>>,
}

impl<T> Component<T> for ParallelCircuit<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        self.components
            .iter()
            .map(|c| c.impedance(freq).finv())
            .sum::<Impedance<T>>()
            .finv()
    }
}

#[derive(Default)]
pub struct SeriesCircuit<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    pub components: Vec<Box<dyn Component<T>>>,
}

impl<T> Component<T> for SeriesCircuit<T>
where
    T: FloatConst + ConstOne + ConstZero,
{
    fn impedance(&self, freq: T) -> Impedance<T> {
        self.components.iter().map(|c| c.impedance(freq)).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::components::rcl::{Capacitor, Resistor};
    use crate::components::Component;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_parallel_r_circuits() {
        let r1 = Resistor { r0: 2.0e3 };
        let r2 = Resistor { r0: 6.0e3 };
        let components: Vec<Box<dyn Component<f64>>> = vec![Box::new(r1), Box::new(r2)];
        let circuit = ParallelCircuit {
            components: components,
        };
        let z = circuit.impedance(100.0);
        assert_approx_eq!(z.re(), 1.5e3);
        assert_approx_eq!(z.im(), 0.0);
    }

    #[test]
    fn test_series_r_circuits() {
        let r1 = Resistor { r0: 2.0e3 };
        let r2 = Resistor { r0: 6.0e3 };
        let components: Vec<Box<dyn Component<f64>>> = vec![Box::new(r1), Box::new(r2)];
        let circuit = SeriesCircuit {
            components: components,
        };
        let z = circuit.impedance(100.0);
        assert_approx_eq!(z.re(), 8.0e3);
        assert_approx_eq!(z.im(), 0.0);
    }

    #[test]
    fn test_parallel_rc_circuits() {
        let r1 = Resistor { r0: 2.0e3 };
        let r2 = Resistor { r0: 6.0e3 };
        let c = Capacitor { c0: 1.0e-5 };
        let components: Vec<Box<dyn Component<f64>>> =
            vec![Box::new(r1), Box::new(r2), Box::new(c)];
        let circuit = ParallelCircuit {
            components: components,
        };
        let z = circuit.impedance(100.0);
        assert_approx_eq!(z.re(), 8.0e3);
        assert_approx_eq!(z.im(), -1.0e3 / std::f64::consts::TAU);
    }

    #[test]
    fn test_series_rc_circuits() {
        let r1 = Resistor { r0: 2.0e3 };
        let r2 = Resistor { r0: 6.0e3 };
        let c = Capacitor { c0: 1.0e-5 };
        let components: Vec<Box<dyn Component<f64>>> =
            vec![Box::new(r1), Box::new(r2), Box::new(c)];
        let circuit = SeriesCircuit {
            components: components,
        };
        let z = circuit.impedance(100.0);
        assert_approx_eq!(z.re(), 8.0e3);
        assert_approx_eq!(z.im(), -1.0e3 / std::f64::consts::TAU);
    }
}