use assert_approx_eq::assert_approx_eq;
use eispec::circuits::{ParallelCircuit, SeriesCircuit};
use eispec::components::cpe::Cpe;
use eispec::components::rlc::Resistor;
use eispec::components::Component;
use eispec::data::ImpedanceDataError;
use eispec::data::{ImpedanceData, ImpedanceDataAccessors};
use eispec::fit::{ImpedanceModel, ModelParameter};
use eispec::newtypes::{Frequency, Impedance};

struct ImpedanceDataWrap(ImpedanceData<f64>);

impl ImpedanceDataAccessors<f64> for ImpedanceDataWrap {
    fn get_freqs(&self) -> &[Frequency<f64>] {
        self.0.get_freqs()
    }
    fn get_zmeas(&self) -> &[Impedance<f64>] {
        self.0.get_zmeas()
    }
    fn get_zerr(&self) -> &[Impedance<f64>] {
        self.0.get_zerr()
    }
    fn get_parameters(&self) -> Option<&[ModelParameter<f64>]> {
        self.0.get_parameters()
    }
    fn set_parameters(&mut self, parameters: Vec<ModelParameter<f64>>) {
        self.0.set_parameters(parameters)
    }
    fn from_csv(filename: &str) -> Result<ImpedanceDataWrap, ImpedanceDataError> {
        ImpedanceData::<f64>::from_csv(filename).map(|a| ImpedanceDataWrap(a))
    }
}

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
        ModelParameter::new(true, Some(0.0), None),
        ModelParameter::new(true, Some(0.0), None),
        ModelParameter::new(true, Some(0.0), Some(1.0)),
        ModelParameter::new(true, Some(10.0), None),
    ]);

    let mut params = [10.0, 1.0e5, 0.1, 10.0];

    let res = data.fit(&mut params).unwrap();

    println!("{:#?}", res);
    println!("{:#?}", params);

    assert_eq!(res.n_par, 4);
    assert_eq!(res.n_free, 4);
    assert_eq!(res.n_func, 15);

    assert_approx_eq!(950.6002232884109, params[0], 0.1);
    assert_approx_eq!(33318.887554247776, params[1], 0.1);
    assert_approx_eq!(0.3498288404822108, params[2], 1e-3);
    assert_approx_eq!(19.740535282535895, params[3], 0.1);
}
