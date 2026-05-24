use assert_approx_eq::assert_approx_eq;
use eispec::circuits::{ParallelCircuit, SeriesCircuit};
use eispec::components::cpe::Cpe;
use eispec::components::rlc::Resistor;
use eispec::components::Component;
use eispec::data::ImpedanceDataError;
use eispec::data::{ImpedanceData, ImpedanceDataAccessors};
use eispec::fit::{ImpedanceModel, ModelParameter, ParameterBounds};
use eispec::impl_impedance_data_accessors;
use eispec::newtypes::{Frequency, Impedance};

#[impl_impedance_data_accessors]
struct ImpedanceDataWrap(ImpedanceData<f64>);

impl ImpedanceModel<f64> for ImpedanceDataWrap {
    fn model(&self, params: &[f64]) -> Box<dyn Component<f64>> {
        let r1 = Resistor::<f64>::new(params[0]);
        let cpe = Cpe::<f64>::new(params[1], params[2]);
        let mut sub = ParallelCircuit::<f64>::new();
        sub.add(r1);
        sub.add(cpe);

        let mut model = SeriesCircuit::<f64>::new();
        let r2 = Resistor::<f64>::new(params[3]);
        model.add(r2);
        model.add(sub);

        Box::new(model)
    }
}

fn main() {
    let mut data = ImpedanceDataWrap::from_csv("examples/simple_example_data.csv").unwrap();

    data.set_parameters(vec![
        ModelParameter::new(10.0, true, ParameterBounds::positive()),
        ModelParameter::new(1.0e5, true, ParameterBounds::positive()),
        ModelParameter::new(0.1, true, ParameterBounds::zero_to_one()),
        ModelParameter::new(20.0, false, ParameterBounds::new(Some(10.0), None)),
    ]);

    let result = data.fit().unwrap();

    println!("{:#?}", data.get_parameters().unwrap());
    println!("{:#?}", result);

    assert_eq!(result.n_par, 4);
    assert_eq!(result.n_free, 3);
    assert_eq!(result.n_func, 15);

    assert_approx_eq!(950.0, result.x[0], 0.2);
    assert_approx_eq!(33333.3, result.x[1], 50.0);
    assert_approx_eq!(0.35, result.x[2], 1e-3);
    assert_eq!(20.0, result.x[3]);
}
