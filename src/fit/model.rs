//!
//! Complex non-linear fitting for impedance data
//!
use num::traits::NumAssign;

use crate::components::Component;
use crate::constants::FloatConst;
use crate::data::ImpedanceDataAccessors;
use crate::fit::config::MPFitConfig;
use crate::fit::enums::{MPFitDone, MPFitError, MPFitInfo};
use crate::fit::mpfit::MPFit;
use crate::fit::status::MPFitStatus;

///
/// Trait to be implemented by the user on some impedance data to be fitted.
///
pub trait ImpedanceModel<T>
where
    T: NumAssign + FloatConst,
    Self: Sized + ImpedanceDataAccessors<T>,
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
            .zip(self.get_freqs().iter())
            .zip(self.get_zmeas().iter())
            .zip(self.get_zerr().iter())
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
    fn fit(&mut self) -> Result<MPFitStatus<T>, MPFitError> {
        let mut init: Vec<T> = self
            .get_parameters()
            .unwrap()
            .iter()
            .map(|p| p.init_value)
            .collect();

        let config = self.config();
        let mut fit = MPFit::new(self, &mut init, &config)?;

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
