//!
//! Jacobian matrix representation used for the complex fitting procedure.
//!
use std::ops::{Index, IndexMut, Range};

use crate::constants::FloatConst;

///
/// Column-major matrix storing the Jacobian.
///
/// Exposes three indexing forms:
/// - `m[k]`       flat (column-major) for tight running-index loops
/// - `m[(i, j)]`  2D element access
/// - `m[a..b]`    flat slice range, for `enorm()` and `iter_mut()`
///
#[derive(Debug, Clone)]
pub struct JacMatrix<T> {
    data: Vec<T>,
    rows: usize,
}

impl<T> JacMatrix<T>
where
    T: FloatConst,
{
    ///
    /// Create and empty matrix.
    ///
    pub fn empty() -> Self {
        Self {
            data: vec![],
            rows: 0,
        }
    }

    ///
    /// Initialize an 0-matrix with the given number of rows and columns.
    ///
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![T::zero(); rows * cols],
            rows,
        }
    }

    ///
    /// Fill all matrix elements with a given value.
    ///
    pub fn fill(&mut self, val: T) {
        self.data.fill(val);
    }

    ///
    /// Iterate through all elements of the column-major matrix.
    ///
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.data.iter()
    }

    ///
    /// Return the jth column of the matrix as a slice.
    ///
    pub fn col(&self, j: usize) -> &[T] {
        let start = j * self.rows;
        &self.data[start..start + self.rows]
    }

    ///
    /// Return the jth column of the matrix as a mutable slice.
    ///
    pub fn col_mut(&mut self, j: usize) -> &mut [T] {
        let start = j * self.rows;
        &mut self.data[start..start + self.rows]
    }

    ///
    /// Swap ath and bth columns of the matrix in place.
    ///
    pub fn swap_cols(&mut self, a: usize, b: usize) {
        if a == b {
            return;
        }
        let (lo, hi) = if a < b { (a, b) } else { (b, a) };
        let (left, right) = self.data.split_at_mut(hi * self.rows);
        let lo_start = lo * self.rows;
        left[lo_start..lo_start + self.rows].swap_with_slice(&mut right[..self.rows]);
    }
}

impl<T> Index<usize> for JacMatrix<T> {
    type Output = T;
    fn index(&self, k: usize) -> &T {
        &self.data[k]
    }
}

impl<T> IndexMut<usize> for JacMatrix<T> {
    fn index_mut(&mut self, k: usize) -> &mut T {
        &mut self.data[k]
    }
}

impl<T> Index<(usize, usize)> for JacMatrix<T> {
    type Output = T;
    fn index(&self, (i, j): (usize, usize)) -> &T {
        &self.data[j * self.rows + i]
    }
}

impl<T> IndexMut<(usize, usize)> for JacMatrix<T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut T {
        &mut self.data[j * self.rows + i]
    }
}

impl<T> Index<Range<usize>> for JacMatrix<T> {
    type Output = [T];
    fn index(&self, r: Range<usize>) -> &[T] {
        &self.data[r]
    }
}

impl<T> IndexMut<Range<usize>> for JacMatrix<T> {
    fn index_mut(&mut self, r: Range<usize>) -> &mut [T] {
        &mut self.data[r]
    }
}

impl<T> IntoIterator for JacMatrix<T> {
    type Item = T;
    type IntoIter = std::vec::IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}
