use assert_approx_eq::assert_approx_eq;
use eispec::circuits::ParallelCircuit;
use eispec::components::dielectric::ColeCole;
use eispec::components::rlc::Capacitor;
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
        let zarc = ColeCole::<f64>::new(params[0], params[1], params[2], params[3]);
        let cp = Capacitor::<f64>::new(params[4]);
        let mut model = ParallelCircuit::<f64>::new();

        model.add(cp);
        model.add(zarc);

        Box::new(model)
    }
}

fn main() {
    let mut data = ImpedanceDataWrap::from_csv("examples/circuit_a_example_data.csv").unwrap();
    data.set_parameters(vec![
        ModelParameter::new(1.2e6, true, ParameterBounds::positive()),
        ModelParameter::new(9.1e5, true, ParameterBounds::positive()),
        ModelParameter::new(0.2, true, ParameterBounds::new(Some(0.0), Some(5.0))),
        ModelParameter::new(0.384, true, ParameterBounds::zero_to_one()),
        ModelParameter::new(1.3e-12, true, ParameterBounds::positive()),
    ]);

    let result = data.fit().unwrap();

    println!("{:#?}", data.get_parameters().unwrap());
    println!("{:#?}", result);

    assert_eq!(result.n_par, 5);
    assert_eq!(result.n_free, 5);
    assert_eq!(result.n_func, 27);

    assert_approx_eq!(996807.8, result.x[0], 0.1);
    assert_approx_eq!(1997064.0, result.x[1], 1.0);
    assert_approx_eq!(0.9863895, result.x[2]);
    assert_approx_eq!(0.29748956, result.x[3]);
    assert_approx_eq!(9.951745e-13, result.x[4], 1e-16);
}
