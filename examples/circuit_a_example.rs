use eispec::circuits::ParallelCircuit;
use eispec::components::dielectric::ColeCole;
use eispec::components::rlc::Capacitor;
use eispec::components::Component;
use eispec::data::ImpedanceDataError;
use eispec::data::{ImpedanceData, ImpedanceDataAccessors};
use eispec::fit::{ImpedanceModel, ModelParameter};
use eispec::impl_impedance_data_accessors;
use eispec::newtypes::{Frequency, Impedance};

#[impl_impedance_data_accessors]
struct ImpedanceDataWrap(ImpedanceData<f64>);

impl ImpedanceModel<f64> for ImpedanceDataWrap {
    fn model(&self, params: &[f64]) -> Box<(dyn Component<f64>)> {
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
        ModelParameter::new(1.2e6, true, Some(0.0), None),
        ModelParameter::new(9.1e5, true, Some(0.0), None),
        ModelParameter::new(0.2, true, Some(0.0), Some(5.0)),
        ModelParameter::new(0.384, true, Some(0.0), Some(1.0)),
        ModelParameter::new(1.3e-12, true, Some(0.0), None),
    ]);

    let res = data.fit().unwrap();

    println!("{:#?}", data.get_parameters().unwrap());
    println!("{:#?}", res);
}
