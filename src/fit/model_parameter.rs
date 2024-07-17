use crate::constants::FloatConst;

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct ModelParameter<T> {
    pub fit: bool,
    pub limit_lower: Option<T>,
    pub limit_upper: Option<T>,
}

impl<T> ModelParameter<T>
where
    T: FloatConst,
{
    pub fn new(fit: bool, limit_lower: Option<T>, limit_upper: Option<T>) -> Self {
        Self {
            fit,
            limit_lower,
            limit_upper,
        }
    }
}
