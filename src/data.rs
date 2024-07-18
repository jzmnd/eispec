use std::io;
use std::num::ParseFloatError;
use std::str::FromStr;
use thiserror::Error;

use crate::constants::FloatConst;
use crate::fit::ModelParameter;
use crate::newtypes::{Frequency, Impedance};

#[derive(Debug, Error)]
pub enum ImpedanceDataError {
    #[error("Unable to parse numeric data from str")]
    StrParseError,
    #[error("Unable to read file")]
    ReadError,
    #[error("Unable read CSV record")]
    CsvError,
}

impl From<io::Error> for ImpedanceDataError {
    fn from(_: io::Error) -> Self {
        ImpedanceDataError::ReadError
    }
}

impl From<csv::Error> for ImpedanceDataError {
    fn from(_: csv::Error) -> Self {
        ImpedanceDataError::CsvError
    }
}

impl From<ParseFloatError> for ImpedanceDataError {
    fn from(_: ParseFloatError) -> Self {
        ImpedanceDataError::StrParseError
    }
}

pub trait ImpedanceDataAccessors<T> {
    ///
    /// Getter method that should return the frequency data.
    ///
    fn get_freqs(&self) -> &[Frequency<T>];
    ///
    /// Getter method that should return the measured impedance data.
    ///
    fn get_zmeas(&self) -> &[Impedance<T>];
    ///
    /// Getter method that should return the experimental error on the
    /// real and imaginary parts of the impedance data.
    ///
    fn get_zerr(&self) -> &[Impedance<T>];
    ///
    /// Getter method that should return the model parameter configs.
    ///
    /// Parameters are expected in the same order as those passed into
    /// the `model(params)` function.
    ///
    fn get_parameters(&self) -> Option<&[ModelParameter<T>]>;
    ///
    /// Set the model parameter configs.
    ///
    fn set_parameters(&mut self, parameters: Vec<ModelParameter<T>>);
    ///
    /// Load the impedance data from a csv file.
    ///
    /// Expected columns: Frequency, Re(Z), Im(Z), Re(Zerr), Im(Zerr)
    ///
    fn from_csv(filename: &str) -> Result<Self, ImpedanceDataError>
    where
        Self: Sized;
}

#[derive(Debug)]
pub struct ImpedanceData<T> {
    pub freqs: Vec<Frequency<T>>,
    pub zmeas: Vec<Impedance<T>>,
    pub zerr: Vec<Impedance<T>>,
    pub parameters: Option<Vec<ModelParameter<T>>>,
}

impl<T> ImpedanceDataAccessors<T> for ImpedanceData<T>
where
    T: FloatConst + FromStr<Err = ParseFloatError>,
{
    fn get_freqs(&self) -> &[Frequency<T>] {
        &self.freqs
    }
    fn get_zmeas(&self) -> &[Impedance<T>] {
        &self.zmeas
    }
    fn get_zerr(&self) -> &[Impedance<T>] {
        &self.zerr
    }
    fn get_parameters(&self) -> Option<&[ModelParameter<T>]> {
        self.parameters.as_deref()
    }
    fn set_parameters(&mut self, parameters: Vec<ModelParameter<T>>) {
        self.parameters = Some(parameters);
    }
    fn from_csv(filename: &str) -> Result<Self, ImpedanceDataError> {
        let mut reader = csv::Reader::from_path(filename)?;
        let mut freqs = Vec::new();
        let mut zmeas = Vec::new();
        let mut zerr = Vec::new();

        for result in reader.records() {
            let record = result?;
            freqs.push(Frequency::new(record[0].parse()?));
            zmeas.push(Impedance::new(record[1].parse()?, record[2].parse()?));
            zerr.push(Impedance::new(record[3].parse()?, record[4].parse()?));
        }
        Ok(ImpedanceData {
            freqs,
            zmeas,
            zerr,
            parameters: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_impedance_data() {
        let data = ImpedanceData::<f64>::from_csv("examples/simple_example_data.csv").unwrap();
        assert_eq!(data.freqs.len(), 15);
        assert_eq!(data.freqs[0].value(), 20.0);
        assert_eq!(data.freqs[14].value(), 1.0e6);
    }

    #[test]
    fn test_set_parameters() {
        let mut data = ImpedanceData::<f64>::from_csv("examples/simple_example_data.csv").unwrap();
        let test_params = vec![ModelParameter::new(true, Some(0.0), None)];
        data.set_parameters(test_params.clone());
        assert_eq!(data.get_parameters().unwrap(), test_params);
    }
}
