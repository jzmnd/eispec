use assert_approx_eq::assert_approx_eq;
use eispec::circuits::{ParallelCircuit, SeriesCircuit};
use eispec::components::cpe::Cpe;
use eispec::components::rlc::Resistor;
use eispec::components::warburg::Warburg;
use eispec::components::Component;
use eispec::data::ImpedanceDataError;
use eispec::data::{ImpedanceData, ImpedanceDataAccessors};
use eispec::fit::config::MPFitConfig;
use eispec::fit::{ImpedanceModel, ModelParameter, ParameterBounds};
use eispec::impl_impedance_data_accessors;
use eispec::newtypes::{Frequency, Impedance};

#[impl_impedance_data_accessors]
struct ImpedanceDataWrap(ImpedanceData<f64>);

impl ImpedanceModel<f64> for ImpedanceDataWrap {
    fn model(&self, params: &[f64]) -> Box<dyn Component<f64>> {
        let r_bulk = Resistor::<f64>::new(params[0]);

        let r_sei = Resistor::<f64>::new(params[1]);
        let cpe_sei = Cpe::<f64>::new(params[2], params[3]);

        let r_ct = Resistor::<f64>::new(params[4]);
        let w = Warburg::<f64>::new(params[5]);
        let cpe_electrode = Cpe::<f64>::new(params[6], params[7]);

        let mut sei = ParallelCircuit::<f64>::new();
        sei.add(r_sei);
        sei.add(cpe_sei);

        let mut ct = ParallelCircuit::<f64>::new();
        let mut s = SeriesCircuit::<f64>::new();
        s.add(r_ct);
        s.add(w);
        ct.add(s);
        ct.add(cpe_electrode);

        let mut model = SeriesCircuit::<f64>::new();
        model.add(r_bulk);
        model.add(sei);
        model.add(ct);

        Box::new(model)
    }

    fn config(&self) -> MPFitConfig<f64> {
        MPFitConfig {
            ftol: 1e-10,
            xtol: 1e-10,
            gtol: 1e-10,
            epsfcn: f64::EPSILON,
            step_factor: 100.0,
            covtol: 1e-14,
            max_iter: 1000,
            max_fev: 0,
            do_user_scale: false,
            no_finite_check: false,
        }
    }
}

fn main() {
    let mut data = ImpedanceDataWrap::from_csv("examples/battery_example_data.csv").unwrap();
    data.set_parameters(vec![
        ModelParameter::new(10.0, true, ParameterBounds::positive()),
        ModelParameter::new(100.0, true, ParameterBounds::positive()),
        ModelParameter::new(1e7, true, ParameterBounds::positive()),
        ModelParameter::new(0.8, true, ParameterBounds::zero_to_one()),
        ModelParameter::new(100.0, true, ParameterBounds::positive()),
        ModelParameter::new(100.0, true, ParameterBounds::positive()),
        ModelParameter::new(1e5, true, ParameterBounds::positive()),
        ModelParameter::new(0.8, true, ParameterBounds::zero_to_one()),
    ]);

    let result = data.fit().unwrap();

    println!("{:#?}", data.get_parameters().unwrap());
    println!("{:#?}", result);

    assert_eq!(result.n_par, 8);
    assert_eq!(result.n_free, 8);
    assert_eq!(result.n_func, 113);

    assert_approx_eq!(30.0, result.x[0], 0.01);
    assert_approx_eq!(75.0, result.x[1], 0.01);
    assert_approx_eq!(25000000.0, result.x[2], 10000.0);
    assert_approx_eq!(0.85, result.x[3], 0.01);
    assert_approx_eq!(150.0, result.x[4], 0.01);
    assert_approx_eq!(200.0, result.x[5], 0.1);
    assert_approx_eq!(333333.3, result.x[6], 1000.0);
    assert_approx_eq!(0.75, result.x[7], 0.01);
}
