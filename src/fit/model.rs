//!
//! Complex non-linear fitting for impedance data
//!
use num::traits::NumAssign;

use crate::components::Component;
use crate::constants::FloatConst;
use crate::fit::config::MPFitConfig;
use crate::fit::enums::{MPFitDone, MPFitError, MPFitInfo};
use crate::fit::mpfit::MPFit;
use crate::fit::status::MPFitStatus;
use crate::fit::ModelParameter;
use crate::newtypes::{Frequency, Impedance};

///
/// Trait to be implemented by the user on some impedance data to be fitted.
///
pub trait ImpedanceModel<T>
where
    T: NumAssign + FloatConst,
    Self: Sized,
{
    ///
    /// Electrochemical model to be fit.
    ///
    /// Should return an electrical component which could be
    /// a single element or multiple elements combined in series
    /// or parallel circuits.
    ///
    fn model(&self, params: &[T]) -> Box<dyn Component<T>>;
    ///
    /// Getter method that should return the frequency data.
    ///
    fn freqs(&self) -> &[Frequency<T>];
    ///
    /// Getter method that should return the measured impedance data.
    ///
    fn zmeas(&self) -> &[Impedance<T>];
    ///
    /// Getter method that should return the experimental error on the
    /// real and imaginary parts of the impedance data.
    ///
    fn zerr(&self) -> &[Impedance<T>];
    ///
    /// Getter method that should return the model parameter configs.
    ///
    /// Parameters are expected in the same order as those passed into
    /// the `model(params)` function.
    ///
    fn parameters(&self) -> Option<&[ModelParameter<T>]>;

    ///
    /// Returns the fit configuration. Can be overidden by the user.
    ///
    fn config(&self) -> MPFitConfig<T> {
        MPFitConfig::default()
    }

    ///
    /// Main evaluation procedure which is called from MPFit.
    ///
    /// The residuals are defined as ```(zmeas[i] - model(freq[i])) / zerr[i]```.
    /// Residuals for the real and imaginary parts of the impedance
    /// are calculated separately and combined into a single value for
    /// the `deviates` slice.
    ///
    fn evaluate(&mut self, params: &[T], deviates: &mut [T]) -> Result<(), MPFitError> {
        let model = self.model(params);

        for (((d, f), zm), ze) in deviates
            .iter_mut()
            .zip(self.freqs().iter())
            .zip(self.zmeas().iter())
            .zip(self.zerr().iter())
        {
            let z = model.impedance(*f);

            let dre = (zm.re() - z.re()) / ze.re();
            let dim = (zm.im() - z.im()) / ze.im();

            *d = (dre.powi(2) + dim.powi(2)).sqrt();
        }
        Ok(())
    }

    ///
    /// Perform the complex non-linear fitting procedure.
    ///
    fn fit(&mut self, init: &mut [T]) -> Result<MPFitStatus<T>, MPFitError> {
        let config = self.config();
        let mut fit = MPFit::new(self, init, &config)?;

        fit.check_config()?;
        fit.parse_parameters()?;
        fit.init_lm()?;

        loop {
            fit.fill_xnew();
            fit.fdjac2()?;
            fit.check_limits();
            fit.qrfac();
            fit.scale();
            fit.transpose();
            if !fit.check_is_finite() {
                return Err(MPFitError::Nan);
            }
            let gnorm = fit.gnorm();
            if gnorm <= config.gtol {
                fit.info = MPFitInfo::ConvergenceDir;
            }
            if fit.info != MPFitInfo::NotDone {
                return fit.terminate();
            }
            if config.max_iter == 0 {
                fit.info = MPFitInfo::MaxIterReached;
                return fit.terminate();
            }
            fit.rescale();
            loop {
                fit.lmpar();
                match fit.iterate(gnorm)? {
                    MPFitDone::Exit => return fit.terminate(),
                    MPFitDone::Inner => continue,
                    MPFitDone::Outer => break,
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuits::{ParallelCircuit, SeriesCircuit};
    use crate::components::cpe::Cpe;
    use crate::components::rlc::Resistor;
    use crate::components::Component;
    use crate::newtypes::{Frequency, Impedance};
    use assert_approx_eq::assert_approx_eq;

    struct ImpedanceData<T> {
        pub freqs: Vec<Frequency<T>>,
        pub zmeas: Vec<Impedance<T>>,
        pub zerr: Vec<Impedance<T>>,
        pub parameters: Option<Vec<ModelParameter<T>>>,
    }

    impl<T> ImpedanceModel<T> for ImpedanceData<T>
    where
        T: NumAssign + FloatConst + Default + 'static,
    {
        fn model(&self, params: &[T]) -> Box<(dyn Component<T> + 'static)> {
            let r1 = Resistor::<T>::new(params[0]);
            let cpe = Cpe::<T>::new(params[1], params[2]);
            let mut sub = ParallelCircuit::<T>::new();
            sub.add(r1);
            sub.add(cpe);

            let mut model = SeriesCircuit::<T>::new();
            let r2 = Resistor::<T>::new(params[3]);
            model.add(r2);
            model.add(sub);

            Box::new(model)
        }

        fn freqs(&self) -> &[Frequency<T>] {
            &self.freqs
        }
        fn zmeas(&self) -> &[Impedance<T>] {
            &self.zmeas
        }
        fn zerr(&self) -> &[Impedance<T>] {
            &self.zerr
        }
        fn parameters(&self) -> Option<&[ModelParameter<T>]> {
            self.parameters.as_deref()
        }
    }

    #[test]
    fn test_fit() {
        let mut data: ImpedanceData<f64> = ImpedanceData {
            freqs: vec![
                Frequency::new(20.0),
                Frequency::new(50.0),
                Frequency::new(100.0),
                Frequency::new(200.0),
                Frequency::new(500.0),
                Frequency::new(1.0e3),
                Frequency::new(2.0e3),
                Frequency::new(5.0e3),
                Frequency::new(1.0e4),
                Frequency::new(2.0e4),
                Frequency::new(5.0e4),
                Frequency::new(1.0e5),
                Frequency::new(2.0e5),
                Frequency::new(5.0e5),
                Frequency::new(1.0e6),
            ],
            zmeas: vec![
                Impedance::new(855.1080956916633, -59.16861452198911),
                Impedance::new(817.2069977881932, -75.30555799831711),
                Impedance::new(781.5684305324787, -87.92476629132373),
                Impedance::new(739.4105089684267, -100.54378026780387),
                Impedance::new(675.0614978734461, -115.72734240629535),
                Impedance::new(618.8644099268437, -125.83212312302642),
                Impedance::new(560.1782231920959, -131.74197352824683),
                Impedance::new(478.4632364664801, -133.32097000356936),
                Impedance::new(416.9898834924106, -130.50558114808908),
                Impedance::new(358.1988032109323, -123.20310907403243),
                Impedance::new(286.9877688694263, -109.80320088682980),
                Impedance::new(240.6860210573763, -98.03711538137328),
                Impedance::new(199.8565324926007, -84.56882738101699),
                Impedance::new(156.6866455660058, -68.77302834015914),
                Impedance::new(129.3131274372006, -57.47773552749856),
            ],
            zerr: vec![Impedance::new(1.0, 1.0); 15],
            parameters: Some(vec![
                ModelParameter::new(true, Some(0.0), None),
                ModelParameter::new(true, Some(0.0), None),
                ModelParameter::new(true, Some(0.0), Some(1.0)),
                ModelParameter::new(true, Some(10.0), None),
            ]),
        };

        let mut init = [10.0, 1.0e5, 0.1, 10.0];

        let res = data.fit(&mut init).unwrap();

        assert_eq!(res.n_par, 4);
        assert_eq!(res.n_free, 4);
        assert_eq!(res.n_func, 15);

        assert_approx_eq!(950.6002232884109, init[0], 0.1);
        assert_approx_eq!(33318.887554247776, init[1], 0.1);
        assert_approx_eq!(0.3498288404822108, init[2], 1e-3);
        assert_approx_eq!(19.740535282535895, init[3], 0.1);
    }
}
