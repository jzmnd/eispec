use crate::constants::FloatConst;

pub fn freq_to_angular<T: FloatConst>(freq: T) -> T {
    T::PI_2 * freq
}
