//! An implementation of 1-dimensional linear interpolation in Rust, similar to MATLAB's `interp1`
//! or NumPy's `numpy.interp`.
//!
//! [`numpy.interp`]: https://numpy.org/doc/stable/reference/generated/numpy.interp.html
//! [`interp1`]: https://www.mathworks.com/help/matlab/ref/double.interp1.html
//!
//! # Example
//!
//! ```
//! use interp::{interp, interp_array, interp_slice, InterpMode};
//!
//! let x = vec![0.0, 0.2, 0.5, 0.8, 1.0];
//! let y = vec![0.0, 1.0, 3.0, 3.5, 4.0];
//!
//! // Interpolate at a single point
//! assert_eq!(interp(&x, &y, 0.35, &InterpMode::default()), 2.0);
//!
//! // Interpolate a slice - alloces a new results Vec<T>
//! let xp = vec![0.1, 0.65, 0.9];
//! assert_eq!(interp_slice(&x, &y, &xp, &InterpMode::default()), vec![0.5, 3.25, 3.75]);
//!
//! // Interpolate an array
//! let xp = [0.1, 0.65, 0.9];
//! assert_eq!(interp_array(&x, &y, &xp, &InterpMode::default()), [0.5, 3.25, 3.75]);
//! ```

#![warn(missing_docs)]
#![allow(unknown_lints)]
#![warn(clippy::all, clippy::pedantic, clippy::cargo)]
#![allow(clippy::many_single_char_names)]

use num_traits::Num;

use itertools::{izip, Itertools};

/// Interpolation method that sets the behaviour of the interpolation functions for points outside
/// the range of sample points.
pub enum InterpMode<T: std::cmp::PartialOrd + Copy> {
    /// Use slope information to linearly extrapolate the value.
    ///
    /// This was the default behaviour in earlier versions of this crate.
    Extrapolate,
    /// Use `y.first()` for `xp <= x.first()`, or `y.last()` for `xp >= x.last()`. This is similar
    /// to the default behaviour or NumPy's [`numpy.interp`] function.
    ///
    /// [`numpy.interp`]: https://numpy.org/doc/stable/reference/generated/numpy.interp.html
    FirstLast,
    /// Use the given constant for values outside the range.
    ///
    /// This is commonly used to return `T::NAN` (assuming `T: Float`) to mirror the MATLAB
    /// [`interp1`] function's behaviour using the `'linear'` method.
    ///
    /// [`interp1`]: https://www.mathworks.com/help/matlab/ref/double.interp1.html
    Constant(T),
}

#[allow(clippy::derivable_impls)] // We need a manual impl for our MSRV
impl<T: std::cmp::PartialOrd + Copy> Default for InterpMode<T> {
    fn default() -> Self {
        InterpMode::Extrapolate
    }
}

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

/// Linearly interpolate the data points given by the `x` and `y` slices at point `xp`,
/// using the interpolation method provided by `mode`.
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
/// use interp::InterpMode;
///
/// let x = vec![0.0, 1.0, 2.0, 3.0];
/// let y = vec![1.0, 3.0, 4.0, 2.0];
///
/// assert_eq!(interp(&x, &y, 1.5, &InterpMode::Extrapolate), 3.5);
/// ```
pub fn interp<T>(x: &[T], y: &[T], xp: T, mode: &InterpMode<T>) -> T
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
        // i = 0 when none found
        let i = prev_index(x, xp).min(min_len - 2);

        let point = m[i] * xp + c[i];
        let x_limits = (&x[0], &x[min_len - 1]);
        let y_limits = (&y[0], &y[min_len - 1]);

        select_outside_point(x_limits, y_limits, &xp, point, mode)
    }
}

/// Linearly interpolate the data points given by the `x` and `y` slices at each of the points in
/// the `xp` slice, using the interpolation method provided by `mode`.
///
/// Returns a `Vec<T>` containing the equivalent y coordinates to each of the x coordinates given
/// by `xp`.
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
/// use interp::InterpMode;
///
/// let x = vec![0.0, 1.0, 2.0, 3.0];
/// let y = vec![1.0, 3.0, 4.0, 2.0];
///
/// let xp = vec![0.5, 2.5, 4.0];
///
/// assert_eq!(interp_slice(&x, &y, &xp, &InterpMode::Extrapolate), vec![2.0, 3.0, 0.0]);
///
/// assert_eq!(interp_slice(&x, &y, &xp, &InterpMode::Constant(1e3)), vec![2.0, 3.0, 1e3]);
///
/// assert_eq!(interp_slice(&x, &y, &xp, &InterpMode::FirstLast), vec![2.0, 3.0, 2.0]);
/// ```
pub fn interp_slice<T>(x: &[T], y: &[T], xp: &[T], mode: &InterpMode<T>) -> Vec<T>
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

        let x_limits = (&x[0], &x[min_len - 1]);
        let y_limits = (&y[0], &y[min_len - 1]);

        xp.iter()
            .map(|&xp| {
                // The index of the x coordinate right before xp. Use min to ensure we don't go out
                // of m's and c's bounds when xp > x.last()
                let i = prev_index(x, xp).min(min_len - 2);

                let point = m[i] * xp + c[i];

                select_outside_point(x_limits, y_limits, &xp, point, mode)
            })
            .collect()
    }
}

/// Linearly interpolate the data points given by the `x` and `y` slices at each of the points in
/// the `xp` slice, using the interpolation method provided by `mode`.
///
/// Returns a `[T; N]` containing the equivalent y coordinates to each of the x coordinates given
/// by `xp`.
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
/// use interp::interp_array;
/// use interp::InterpMode;
///
/// let x = [0.0, 1.0, 2.0, 3.0];
/// let y = [1.0, 3.0, 4.0, 2.0];
///
/// let xp = [0.5, 2.5, 4.0];
///
/// assert_eq!(interp_array(&x, &y, &xp, &InterpMode::Extrapolate), [2.0, 3.0, 0.0]);
/// ```
pub fn interp_array<T, const N: usize>(
    x: &[T],
    y: &[T],
    xp: &[T; N],
    mode: &InterpMode<T>,
) -> [T; N]
where
    T: Num + PartialOrd + Copy,
{
    // The min-length of the x and y vectors. We ignore additional entries in either vec.
    let min_len = std::cmp::min(x.len(), y.len());

    if min_len == 0 {
        [T::zero(); N]
    } else if min_len == 1 {
        [y[0]; N]
    } else {
        // Difference between subsequent x and y coordinate values
        let dx = deltas(&x[..min_len]);
        let dy = deltas(&y[..min_len]);

        // Slope between subsequent points
        let m = slopes(&dx, &dy);

        // Intercept of the line between adjacent points
        let c = intercepts(x, y, &m);

        let x_limits = (&x[0], &x[min_len - 1]);
        let y_limits = (&y[0], &y[min_len - 1]);

        xp.map(|xp| {
            // The index of the x coordinate right before xp. Use min to ensure we don't go out
            // of m's and c's bounds when xp > x.last()
            let i = prev_index(x, xp).min(min_len - 2);

            let point = m[i] * xp + c[i];

            select_outside_point(x_limits, y_limits, &xp, point, mode)
        })
    }
}

fn select_outside_point<T>(
    x_limits: (&T, &T),
    y_limits: (&T, &T),
    xp: &T,
    default: T,
    mode: &InterpMode<T>,
) -> T
where
    T: Num + PartialOrd + Copy,
{
    if xp < x_limits.0 {
        match mode {
            InterpMode::Extrapolate => default,
            InterpMode::FirstLast => *y_limits.0,
            InterpMode::Constant(val) => *val,
        }
    } else if xp > x_limits.1 {
        match mode {
            InterpMode::Extrapolate => default,
            InterpMode::FirstLast => *y_limits.1,
            InterpMode::Constant(val) => *val,
        }
    } else {
        default
    }
}

#[cfg(test)]
mod tests {
    use isclose::assert_is_close;
    use num_traits::float::FloatCore;

    use super::*;

    fn vec_compare<T>(v1: &[T], v2: &[T]) -> bool
    where
        T: Num + PartialOrd + Copy + FloatCore,
    {
        (v1.len() == v2.len()) && izip!(v1, v2).all(|(a, b)| (a.is_nan() && b.is_nan()) || (a == b))
    }

    #[test]
    fn test_deltas() {
        let p = vec![0.0, 1.0, 3.5, 2.5, 6.0];
        let delta = vec![1.0, 2.5, -1.0, 3.5];

        let result = deltas(&p);

        assert_eq!(result.len(), delta.len());

        for (r, d) in izip!(result, delta) {
            assert_is_close!(r, d);
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
            assert_is_close!(r, m);
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
            assert_is_close!(r, c);
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
        assert_is_close!(interp(&[], &[], 2.0, &InterpMode::Extrapolate), 0.0);

        assert_is_close!(interp(&[1.0], &[2.0], 2.0, &InterpMode::Extrapolate), 2.0);

        let x = vec![0.0, 1.0, 2.0, 3.0, 4.5];
        let y = vec![0.0, 2.0, 5.0, 3.0, 2.0];

        assert_is_close!(interp(&x, &y, 2.5, &InterpMode::Extrapolate), 4.0);
        assert_is_close!(interp(&x, &y, -1.0, &InterpMode::Extrapolate), -2.0);
        assert_is_close!(interp(&x, &y, 7.5, &InterpMode::Extrapolate), 0.0);
    }

    #[test]
    fn test_interp_slice() {
        assert_eq!(
            interp_slice(&[], &[], &[2.0], &InterpMode::Extrapolate),
            vec![0.0]
        );

        assert_eq!(
            interp_slice(&[1.0], &[2.0], &[2.0], &InterpMode::Extrapolate),
            vec![2.0]
        );

        let x = vec![0.0, 1.0, 2.0, 3.0, 4.5];
        let y = vec![0.0, 2.0, 5.0, 3.0, 2.0];

        assert_eq!(interp_slice(&x, &y, &[], &InterpMode::Extrapolate), vec![]);

        let xp = vec![2.5, -1.0, 7.5];
        let result = vec![4.0, -2.0, 0.0];

        assert_eq!(interp_slice(&x, &y, &xp, &InterpMode::Extrapolate), result);
    }

    #[test]
    fn test_interp_slice_constant_outside_bounds() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![1.0, 3.0, 4.0, 2.0];

        assert_eq!(
            interp_slice(&x, &y, &[], &InterpMode::Constant(f64::NAN)),
            vec![]
        );

        let xp = vec![0.5, 2.5, 4.0];

        let expected = vec![2.0, 3.0, 1e3];
        assert_eq!(
            interp_slice(&x, &y, &xp, &InterpMode::Constant(1e3)),
            expected
        );

        let x = vec![0.0, 1.0, 2.0, 3.0, 4.5];
        let y = vec![0.0, 2.0, 5.0, 3.0, 2.0];

        let xp = vec![2.5, -1.0, 7.5];
        let result = interp_slice(&x, &y, &xp, &InterpMode::Constant(f64::NAN));
        let expected = vec![4.0, f64::NAN, f64::NAN];

        assert!(vec_compare(&result, &expected));
    }

    #[test]
    fn test_interp_slice_first_last() {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.5];
        let y = vec![0.0, 2.0, 5.0, 3.0, 2.0];

        assert_eq!(interp_slice(&x, &y, &[], &InterpMode::FirstLast), vec![]);

        let xp = vec![2.5, -1.0, 7.5];
        let result = vec![4.0, 0.0, 2.0];

        assert_eq!(interp_slice(&x, &y, &xp, &InterpMode::FirstLast), result);
    }

    #[test]
    fn test_interp_array() {
        let result = interp_array(&[], &[], &[2.0], &InterpMode::Extrapolate);
        for (act, exp) in izip!(result, [0.0]) {
            assert_is_close!(act, exp);
        }

        let result = interp_array(&[1.0], &[2.0], &[2.0], &InterpMode::Extrapolate);
        for (act, exp) in izip!(result, [2.0]) {
            assert_is_close!(act, exp);
        }

        let x = vec![0.0, 1.0, 2.0, 3.0, 4.5];
        let y = vec![0.0, 2.0, 5.0, 3.0, 2.0];
        let result = interp_array(&x, &y, &[], &InterpMode::Extrapolate);
        for (act, exp) in izip!(result, [] as [f64; 0]) {
            assert_is_close!(act, exp);
        }

        let xp = [2.5, -1.0, 7.5];
        let expected = [4.0, -2.0, 0.0];
        let result = interp_array(&x, &y, &xp, &InterpMode::Extrapolate);
        for (act, exp) in izip!(result, expected) {
            assert_is_close!(act, exp);
        }
    }
}
