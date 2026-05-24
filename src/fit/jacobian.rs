use crate::constants::FloatConst;
use std::ops::{Index, IndexMut, Range};

/// Column-major matrix storing the Jacobian (and in-place factorisations of it).
///
/// Exposes three indexing forms:
/// - `m[k]`         — flat (column-major) for tight running-index loops
/// - `m[(i, j)]`    — 2D element access
/// - `m[a..b]`      — flat slice range, for `enorm()` and `iter_mut()`
#[derive(Debug, Clone)]
pub struct JacMatrix<T> {
    data: Vec<T>,
    rows: usize,
}

impl<T> JacMatrix<T>
where
    T: FloatConst
{
    pub fn empty() -> Self {
        Self {
            data: vec![],
            rows: 0,
        }
    }

    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![T::zero(); rows * cols],
            rows,
        }
    }

    pub fn fill(&mut self, val: T) {
        self.data.fill(val);
    }

    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.data.iter()
    }

    pub fn col(&self, j: usize) -> &[T] {
        let start = j * self.rows;
        &self.data[start..start + self.rows]
    }

    pub fn col_mut(&mut self, j: usize) -> &mut [T] {
        let start = j * self.rows;
        &mut self.data[start..start + self.rows]
    }

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
