//!
//! Core circuit components to combine other components in series or parallel
//!
use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::{Frequency, Impedance};

#[derive(Default)]
pub struct ParallelCircuit<T> {
    pub components: Vec<Box<dyn Component<T>>>,
}

impl<T> ParallelCircuit<T>
where
    T: FloatConst + Default,
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, component: impl Component<T> + 'static) {
        self.components.push(Box::new(component));
    }
}

impl<T> Component<T> for ParallelCircuit<T>
where
    T: FloatConst,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        self.components
            .iter()
            .map(|c| c.impedance(freq).finv())
            .sum::<Impedance<T>>()
            .finv()
    }
}

#[derive(Default)]
pub struct SeriesCircuit<T> {
    pub components: Vec<Box<dyn Component<T>>>,
}

impl<T> SeriesCircuit<T>
where
    T: FloatConst + Default,
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, component: impl Component<T> + 'static) {
        self.components.push(Box::new(component));
    }
}

impl<T> Component<T> for SeriesCircuit<T>
where
    T: FloatConst,
{
    fn impedance(&self, freq: Frequency<T>) -> Impedance<T> {
        self.components.iter().map(|c| c.impedance(freq)).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::components::rlc::{Capacitor, Resistor};
    use crate::components::Component;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_parallel_r_circuits() {
        let r1 = Resistor::new(2.0e3);
        let r2 = Resistor::new(6.0e3);
        let f = Frequency::new(100.0);
        let mut circuit = ParallelCircuit::<f64>::new();
        circuit.add(r1);
        circuit.add(r2);
        let z = circuit.impedance(f);
        assert_approx_eq!(z.re(), 1.5e3);
        assert_approx_eq!(z.im(), 0.0);
    }

    #[test]
    fn test_series_r_circuits() {
        let r1 = Resistor::new(2.0e3);
        let r2 = Resistor::new(6.0e3);
        let f = Frequency::new(100.0);
        let mut circuit = SeriesCircuit::<f64>::new();
        circuit.add(r1);
        circuit.add(r2);
        let z = circuit.impedance(f);
        assert_approx_eq!(z.re(), 8.0e3);
        assert_approx_eq!(z.im(), 0.0);
    }

    #[test]
    fn test_parallel_rc_circuits() {
        let r1 = Resistor::new(2.0e3);
        let r2 = Resistor::new(6.0e3);
        let c = Capacitor::new(1.0e-5);
        let f = Frequency::new(100.0);
        let mut circuit = ParallelCircuit::<f64>::new();
        circuit.add(r1);
        circuit.add(r2);
        circuit.add(c);
        let z = circuit.impedance(f);
        assert_approx_eq!(z.re(), 1.669886958e1);
        assert_approx_eq!(z.im(), -1.573831380e2);
    }

    #[test]
    fn test_series_rc_circuits() {
        let r1 = Resistor::new(2.0e3);
        let r2 = Resistor::new(6.0e3);
        let c = Capacitor::new(1.0e-5);
        let f = Frequency::new(100.0);
        let mut circuit = SeriesCircuit::<f64>::new();
        circuit.add(r1);
        circuit.add(r2);
        circuit.add(c);
        let z = circuit.impedance(f);
        assert_approx_eq!(z.re(), 8.0e3);
        assert_approx_eq!(z.im(), -1.591549431e2);
    }

    #[test]
    fn test_r_rc_circuits() {
        let r1 = Resistor::new(2.0e3);
        let r2 = Resistor::new(6.0e3);
        let c = Capacitor::new(1.0e-5);
        let f = Frequency::new(100.0);
        let mut subcircuit = ParallelCircuit::<f64>::new();
        subcircuit.add(r2);
        subcircuit.add(c);
        let mut circuit = SeriesCircuit::<f64>::new();
        circuit.add(r1);
        circuit.add(subcircuit);
        let z = circuit.impedance(f);
        assert_approx_eq!(z.re(), 2.004218748e3);
        assert_approx_eq!(z.im(), -1.590430373e2);
    }
}
