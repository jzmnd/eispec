//!
//! Copied and modified from rmpfit crate Copyright (c) Vadim Dyadkin
//! Rust implementation of [CMPFIT](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html)
//!
use num::traits::NumAssign;

use crate::constants::FloatConst;

///
/// Given an n-vector x, this function calculates the euclidean norm of x.
///
/// The euclidean norm is computed by accumulating the sum of
/// squares in three different sums. The sums of squares for the
/// small and large components are scaled so that no overflows
/// occur. non-destructive underflows are permitted. Underflows
/// and overflows do not occur in the computation of the unscaled
/// sum of squares for the intermediate components.
/// The definitions of small, intermediate and large components
/// depend on two constants, RDWARF and RGIANT. The main
/// restrictions on these constants are that RDWARF**2 not
/// underflow and RGIANT**2 not overflow. The constants
/// given here are suitable for every known computer.
///
pub trait ENorm<T> {
    fn enorm(&self) -> T;
}

impl<T> ENorm<T> for [T]
where
    T: NumAssign + FloatConst,
{
    fn enorm(&self) -> T {
        let mut s1 = T::zero();
        let mut s2 = T::zero();
        let mut s3 = T::zero();
        let mut x1max = T::zero();
        let mut x3max = T::zero();
        let agiant = T::MP_RGIANT / T::from(self.len()).unwrap();

        for &val in self {
            let xabs = val.abs();
            if xabs > T::MP_RDWARF && xabs < agiant {
                // Sum for intermediate components.
                s2 += xabs.powi(2);
            } else if xabs > T::MP_RDWARF {
                // Sum for large components.
                if xabs > x1max {
                    s1 = T::one() + s1 * (x1max / xabs).powi(2);
                    x1max = xabs;
                } else {
                    s1 += (xabs / x1max).powi(2);
                }
            } else if xabs > x3max {
                // Sum for small components.
                s3 = T::one() + s3 * (x3max / xabs).powi(2);
                x3max = xabs;
            } else if xabs != T::zero() {
                s3 += (xabs / x3max).powi(2);
            }
        }
        // Calculation of norm.
        if s1 != T::zero() {
            x1max * (s1 + (s2 / x1max) / x1max).sqrt()
        } else if s2 != T::zero() {
            if s2 >= x3max {
                s2 * (T::one() + (x3max / s2) * (x3max * s3))
            } else {
                x3max * ((s2 / x3max) + (x3max * s3))
            }
            .sqrt()
        } else {
            x3max * s3.sqrt()
        }
    }
}
