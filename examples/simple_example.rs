use assert_approx_eq::assert_approx_eq;
use eispec::circuits::{ParallelCircuit, SeriesCircuit};
use eispec::components::cpe::Cpe;
use eispec::components::rlc::Resistor;
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
        ModelParameter::new(10.0, true, Some(0.0), None),
        ModelParameter::new(1.0e5, true, Some(0.0), None),
        ModelParameter::new(0.1, true, Some(0.0), Some(1.0)),
        ModelParameter::new(20.0, false, Some(10.0), None),
    ]);

    let res = data.fit().unwrap();

    println!("{:#?}", data.get_parameters().unwrap());
    println!("{:#?}", res);

    assert_eq!(res.n_par, 4);
    assert_eq!(res.n_free, 3);
    assert_eq!(res.n_func, 15);

    assert_approx_eq!(950.0, res.x[0], 0.2);
    assert_approx_eq!(33333.3, res.x[1], 50.0);
    assert_approx_eq!(0.35, res.x[2], 1e-3);
    assert_eq!(20.0, res.x[3]);
}
