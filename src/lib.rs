//! A Rust implementation of Matlab's `interp1` function for linear interpolation.
//!
//! # Example
//!
//! ```
//! use interp::interp;
//!
//! let x = vec![0.0, 0.2, 0.5, 0.8, 1.0];
//! let y = vec![0.0, 1.0, 3.0, 3.5, 4.0];
//!
//! assert_eq!(interp(&x, &y, 0.35), 2.0);
//! ```

#![warn(missing_docs)]
#![allow(unknown_lints)]
#![warn(clippy::all, clippy::pedantic, clippy::cargo)]
#![allow(clippy::many_single_char_names)]

use num::Num;

use itertools::{izip, Itertools};

/// Finds the delta between adjacent entries in the slice `p`.
///
/// Returns a [`Vec<T>`] which is one element shorter than `p.len()` containing the difference
/// between each pair of values.
///
/// If `p.len()` is one or less, returns an empty vector.
///
/// # Example
///
/// ```ignore
/// let p = vec![0.0, 1.0, 3.0, 3.5];
///
/// assert_eq!(deltas(&p), vec![1.0, 2.0, 0.5]);
/// ```
#[inline]
fn deltas<T>(p: &[T]) -> Vec<T>
where
    T: Num + Copy,
{
    p.iter().tuple_windows().map(|(&p1, &p2)| p2 - p1).collect()
}

/// Finds the slope of a line segment given two slices containing the x and y deltas.
///
/// Parameters `dx` and `dy` are the differences between adjacent x and y coordinates, respectively.
/// Returns a [`Vec<T>`] containing the slope for each segment.
///
/// If the lengths of `dx` and `dy` are not equal, only the number of elements in the shorter slice
/// are considered; excess elements are ignored.
///
/// # Example
///
/// ```ignore
/// let dx = vec![1.0, 2.0, 0.5];
/// let dy = vec![1.0; 3];
///
/// assert_eq!(slopes(&dx, &dy), vec![1.0, 0.5, 2.0]);
/// ```
#[inline]
fn slopes<T>(dx: &[T], dy: &[T]) -> Vec<T>
where
    T: Num + Copy,
{
    izip!(dx, dy).map(|(&dx, &dy)| dy / dx).collect()
}

/// Finds the y-intercept of line segments given by the points `x` and `y`, and slopes `m`.
///
/// Returns a `Vec<T>` of the same length containing each of the intercepts.
///
/// If the lengths of `x`, `y`, and `m` are not equal, only the number of elements in the shortest
/// slice are considered; excess elements are ignored.
///
/// # Example
///
/// ```ignore
/// let x = vec![0.0, 1.0, 3.5];
/// let y = vec![0.0, 1.0, 6.0];
/// let slope = vec![1.0, 2.0, -0.5];
///
/// assert_eq!(intercepts(&x, &y, &slope), vec![0.0, -1.0, 7.75])
/// ```
#[inline]
fn intercepts<T>(x: &[T], y: &[T], m: &[T]) -> Vec<T>
where
    T: Num + Copy,
{
    izip!(x, y, m).map(|(&x, &y, &m)| y - x * m).collect()
}

/// Finds the index of the value in `x` just before `xp`.
///
/// If the values in `x` are not strictly increasing, the first possible result is returned.
///
/// Returns 0 if there are no elements in `x` which are less than `xp`.
///
/// # Example
///
/// ```ignore
/// let x = vec![0.0, 1.0, 3.0, 4.5];
///
/// assert_eq!(prev_index(&x, 3.5), 2);
/// ```
#[inline]
fn prev_index<T>(x: &[T], xp: T) -> usize
where
    T: Num + PartialOrd + Copy,
{
    x.iter()
        .take_while(|&&x| x < xp)
        .enumerate()
        .last()
        .map_or(0, |(i, _)| i)
}

/// Linearly interpolate the data points given by the `x` and `y` slices at point `xp`.
///
/// Returns the equivalent y coordinate to the x coordinate given by `xp`.
///
/// If the lengths of `x` and `y` differ, only the number of elements in the shorter slice are
/// considered; excess elements are ignored.
///
/// If the length of either `x` or `y` is 0, 0 is returned. If the length of either is 1, `y[0]` is
/// returned. If both are 2 elements or longer the interpolation is performed as expected.
///
/// # Example
///
/// ```
/// use interp::interp;
///
/// let x = vec![0.0, 1.0, 2.0, 3.0];
/// let y = vec![1.0, 3.0, 4.0, 2.0];
///
/// assert_eq!(interp(&x, &y, 1.5), 3.5);
/// ```
pub fn interp<T>(x: &[T], y: &[T], xp: T) -> T
where
    T: Num + PartialOrd + Copy,
{
    // The min-length of the x and y vectors. We ignore additional entries in either vec.
    let min_len = std::cmp::min(x.len(), y.len());

    if min_len == 0 {
        T::zero()
    } else if min_len == 1 {
        y[0]
    } else {
        // Difference between subsequent x and y coordinate values
        let dx = deltas(&x[..min_len]);
        let dy = deltas(&y[..min_len]);

        // Slope between subsequent points
        let m = slopes(&dx, &dy);

        // Intercept of the line between adjacent points
        let c = intercepts(x, y, &m);

        // The index of the x coordinate right before xp
        let i = prev_index(x, xp).min(min_len - 2);

        m[i] * xp + c[i]
    }
}

/// Linearly interpolate the data points given by the `x` and `y` slices at each of the points in
/// the `xp` slice.
///
/// Returns the equivalent y coordinates to each of the x coordinates given by `xp`.
///
/// This is equivalent to running [`interp`] iteratively for each value in `xp`, but more efficient
/// as intermediate calculations are not repeated.
///
/// If the lengths of `x` and `y` differ, only the number of elements in the shorter slice are
/// considered; excess elements are ignored.
///
/// If the length of either `x` or `y` is 0, 0 is returned for each `xp`. If the length of either
/// is 1, `y[0]` is returned. If both are 2 elements or longer the interpolations are performed as
/// expected.
///
/// # Example
///
/// ```
/// use interp::interp_slice;
///
/// let x = vec![0.0, 1.0, 2.0, 3.0];
/// let y = vec![1.0, 3.0, 4.0, 2.0];
///
/// let xp = vec![0.5, 2.5, 4.0];
///
/// assert_eq!(interp_slice(&x, &y, &xp), vec![2.0, 3.0, 0.0]);
/// ```
pub fn interp_slice<T>(x: &[T], y: &[T], xp: &[T]) -> Vec<T>
where
    T: Num + PartialOrd + Copy,
{
    // The min-length of the x and y vectors. We ignore additional entries in either vec.
    let min_len = std::cmp::min(x.len(), y.len());

    if min_len == 0 {
        vec![T::zero(); xp.len()]
    } else if min_len == 1 {
        vec![y[0]; xp.len()]
    } else {
        // Difference between subsequent x and y coordinate values
        let dx = deltas(&x[..min_len]);
        let dy = deltas(&y[..min_len]);

        // Slope between subsequent points
        let m = slopes(&dx, &dy);

        // Intercept of the line between adjacent points
        let c = intercepts(x, y, &m);

        xp.iter()
            .map(|&xp| {
                // The index of the x coordinate right before xp. Use min to ensure we don't go out
                // of m's and c's bounds when xp > x.last()
                let i = prev_index(x, xp).min(min_len - 2);

                m[i] * xp + c[i]
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deltas() {
        let p = vec![0.0, 1.0, 3.5, 2.5, 6.0];
        let delta = vec![1.0, 2.5, -1.0, 3.5];

        let result = deltas(&p);

        assert_eq!(result.len(), delta.len());

        for (r, d) in izip!(result, delta) {
            assert_eq!(r, d);
        }
    }

    #[test]
    fn test_slopes() {
        let dx = vec![1.0, 2.5, 2.0, 1.0];
        let dy = vec![1.0, 5.0, -1.0, 3.5];
        let slope = vec![1.0, 2.0, -0.5, 3.5];

        let result = slopes(&dx, &dy);

        assert_eq!(result.len(), slope.len());

        for (r, m) in izip!(result, slope) {
            assert_eq!(r, m);
        }
    }

    #[test]
    fn test_intercepts() {
        let x = vec![0.0, 1.0, 3.5, 2.5, 6.0];
        let y = vec![0.0, 1.0, 6.0, 5.0, 8.5];
        let slope = vec![1.0, 2.0, -0.5, 3.5, 1.0];
        let intercept = vec![0.0, -1.0, 7.75, -3.75, 2.5];

        let result = intercepts(&x, &y, &slope);

        assert_eq!(result.len(), intercept.len());

        for (r, c) in izip!(result, intercept) {
            assert_eq!(r, c);
        }
    }

    #[test]
    fn test_prev_index() {
        let x = vec![0.0, 1.0, 3.5, 2.5, 6.0];

        let result = prev_index(&x, 3.0);
        assert_eq!(result, 1);

        let result = prev_index(&x, 4.0);
        assert_eq!(result, 3);
    }

    #[test]
    fn test_interp() {
        assert_eq!(interp(&[], &[], 2.0), 0.0);

        assert_eq!(interp(&[1.0], &[2.0], 2.0), 2.0);

        let x = vec![0.0, 1.0, 2.0, 3.0, 4.5];
        let y = vec![0.0, 2.0, 5.0, 3.0, 2.0];

        assert_eq!(interp(&x, &y, 2.5), 4.0);
        assert_eq!(interp(&x, &y, -1.0), -2.0);
        assert_eq!(interp(&x, &y, 7.5), 0.0);
    }

    #[test]
    fn test_interp_slice() {
        assert_eq!(interp_slice(&[], &[], &[2.0]), vec![0.0]);

        assert_eq!(interp_slice(&[1.0], &[2.0], &[2.0]), vec![2.0]);

        let x = vec![0.0, 1.0, 2.0, 3.0, 4.5];
        let y = vec![0.0, 2.0, 5.0, 3.0, 2.0];

        assert_eq!(interp_slice(&x, &y, &[]), vec![]);

        let xp = vec![2.5, -1.0, 7.5];
        let result = vec![4.0, -2.0, 0.0];

        assert_eq!(interp_slice(&x, &y, &xp), result);
    }
}
