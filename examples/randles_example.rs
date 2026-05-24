use assert_approx_eq::assert_approx_eq;
use eispec::circuits::{ParallelCircuit, SeriesCircuit};
use eispec::components::rlc::{Capacitor, Resistor};
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
        let rct = Resistor::<f64>::new(params[0]);
        let cdl = Capacitor::<f64>::new(params[1]);
        let mut sub = ParallelCircuit::<f64>::new();
        sub.add(rct);
        sub.add(cdl);

        let mut model = SeriesCircuit::<f64>::new();
        let rs = Resistor::<f64>::new(params[2]);
        model.add(rs);
        model.add(sub);

        Box::new(model)
    }
}

fn main() {
    let mut data = ImpedanceDataWrap::from_csv("examples/randles_example_data.csv").unwrap();
    data.set_parameters(vec![
        ModelParameter::new(100.0, true, ParameterBounds::positive()),
        ModelParameter::new(2e-6, true, ParameterBounds::positive()),
        ModelParameter::new(500.0, true, ParameterBounds::positive()),
    ]);

    let result = data.fit().unwrap();

    println!("{:#?}", data.get_parameters().unwrap());
    println!("{:#?}", result);

    assert_eq!(result.n_par, 3);
    assert_eq!(result.n_free, 3);
    assert_eq!(result.n_func, 81);

    assert_approx_eq!(200.0, result.x[0], 0.01);
    assert_approx_eq!(1e-6, result.x[1]);
    assert_approx_eq!(1000.0, result.x[2], 0.01);
}
