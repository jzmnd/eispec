use crate::constants::FloatConst;

///
/// Definition of lower/upper bounds on a parameter.
///
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ParameterBounds<T> {
    pub lower: Option<T>,
    pub upper: Option<T>,
}

impl<T> Default for ParameterBounds<T> {
    fn default() -> Self {
        Self {
            lower: None,
            upper: None,
        }
    }
}

impl<T> ParameterBounds<T>
where
    T: FloatConst,
{
    pub fn new(lower: Option<T>, upper: Option<T>) -> Self {
        Self { lower, upper }
    }

    pub fn positive() -> Self {
        Self {
            lower: Some(T::zero()),
            upper: None,
        }
    }

    pub fn negative() -> Self {
        Self {
            lower: None,
            upper: Some(T::zero()),
        }
    }

    pub fn zero_to_one() -> Self {
        Self {
            lower: Some(T::zero()),
            upper: Some(T::one()),
        }
    }
}

///
/// Model parameter configurations.
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct ModelParameter<T> {
    /// Initial parameter value
    pub init_value: T,
    /// Whether to fit the parameter or hold fixed
    pub fit: bool,
    /// Lower/upper limits on parameter, unbounded if None
    pub bounds: ParameterBounds<T>,
}

impl<T> ModelParameter<T>
where
    T: FloatConst,
{
    pub fn new(init_value: T, fit: bool, bounds: ParameterBounds<T>) -> Self {
        Self {
            init_value,
            fit,
            bounds,
        }
    }
}
