use crate::constants::FloatConst;

///
/// Model parameter configurations.
///
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct ModelParameter<T> {
    /// Initial parameter value
    pub init_value: T,
    /// Whether to fit the parameter or hold fixed
    pub fit: bool,
    /// Lower limit on parameter, unbounded if None
    pub limit_lower: Option<T>,
    /// Upper limit on parameter, unbounded if None
    pub limit_upper: Option<T>,
}

impl<T> ModelParameter<T>
where
    T: FloatConst,
{
    pub fn new(init_value: T, fit: bool, limit_lower: Option<T>, limit_upper: Option<T>) -> Self {
        Self {
            init_value,
            fit,
            limit_lower,
            limit_upper,
        }
    }
}
