//!
//! Complex non-linear fitting for impedance data
//!
use num::traits::NumAssign;

use crate::components::Component;
use crate::constants::FloatConst;
use crate::data::ImpedanceDataAccessors;
use crate::fit::config::MPFitConfig;
use crate::fit::enums::MPFitError;
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
    /// Main evaluation procedure which is called from `MPFit`.
    ///
    /// Writes `2 * n` weighted residuals into the `residuals` slice, where
    /// `n` is the number of frequency points. The layout is stacked: the
    /// first `n` entries are real-part residuals `(zmeas.re - z.re) / zerr.re`,
    /// the next `n` are the imaginary-part residuals. Keeping the parts
    /// separate (and signed) preserves the gradient direction.
    ///
    fn evaluate(&mut self, params: &[T], residuals: &mut [T]) -> Result<(), MPFitError> {
        let model = self.model(params);
        let n = self.get_freqs().len();
        let (re_residuals, im_residuals) = residuals.split_at_mut(n);

        for ((((r_re, r_im), f), zm), ze) in re_residuals
            .iter_mut()
            .zip(im_residuals.iter_mut())
            .zip(self.get_freqs().iter())
            .zip(self.get_zmeas().iter())
            .zip(self.get_zerr().iter())
        {
            let z = model.impedance(*f);
            *r_re = (zm.re() - z.re()) / ze.re();
            *r_im = (zm.im() - z.im()) / ze.im();
        }
        Ok(())
    }

    ///
    /// Perform the complex non-linear fitting procedure.
    ///
    fn fit(&mut self) -> Result<MPFitStatus<T>, MPFitError> {
        let init: Vec<T> = self
            .get_parameters()
            .unwrap()
            .iter()
            .map(|p| p.init_value)
            .collect();

        let config = self.config();
        let mut fit = MPFit::try_new(self, init, &config)?;

        let nfree = fit.nfree();
        let mut rdiag = vec![T::zero(); nfree];
        let mut acnorm = vec![T::zero(); nfree];
        let mut step = vec![T::zero(); nfree];

        loop {
            fit.fill_xnew();
            fit.fdjac2()?;
            fit.check_limits();
            fit.qrfac(&mut rdiag, &mut acnorm);
            fit.scale(&acnorm);
            fit.transpose(&rdiag);
            if !fit.check_is_finite() {
                return Err(MPFitError::Nan);
            }
            fit.check_convergence_ortho(&acnorm);
            fit.check_no_iter();
            if fit.is_done() {
                return fit.terminate();
            }
            let gnorm = fit.gnorm(&acnorm);
            fit.rescale(&acnorm);
            loop {
                fit.lmpar(&mut step);
                let accepted = fit.iterate(gnorm, &mut step)?;
                if fit.is_done() {
                    return fit.terminate();
                }
                if accepted {
                    break;
                }
            }
        }
    }
}
