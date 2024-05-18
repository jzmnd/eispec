use num::traits::{ConstOne, ConstZero};

use crate::components::Component;
use crate::constants::FloatConst;
use crate::newtypes::Impedance;

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
