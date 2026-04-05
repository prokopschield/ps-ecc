use crate::Polynomial;

impl Polynomial {
    /// Evaluates a polynomial from multiple coefficient slices at the given point.
    ///
    /// Slices are concatenated in order, with the first slice containing the lowest degree terms.
    #[must_use]
    pub fn eval_coefficient_slices_at(slices: &[&[u8]], x: u8) -> u8 {
        let iter = slices.iter().flat_map(|s| s.iter()).copied();

        Self::eval_coefficient_iter_at(iter, x)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        finite_field::{add, mul},
        Polynomial,
    };

    /// Computes expected value using naive polynomial evaluation.
    fn naive_eval(coefficients: impl Iterator<Item = u8>, x: u8) -> u8 {
        let mut result = 0u8;
        let mut x_pow = 1u8;

        for c in coefficients {
            result = add(result, mul(c, x_pow));
            x_pow = mul(x_pow, x);
        }

        result
    }

    /// Asserts that `eval_coefficient_slices_at` matches naive evaluation for all x in 0..=255.
    fn assert_exhaustive(slices: &[&[u8]]) {
        for x in 0..=255u8 {
            let iter = slices.iter().flat_map(|s| s.iter()).copied();

            assert_eq!(
                Polynomial::eval_coefficient_slices_at(slices, x),
                naive_eval(iter, x)
            );
        }
    }

    #[test]
    fn empty_slice_list() {
        assert_exhaustive(&[]);
    }

    #[test]
    fn single_empty_slice() {
        assert_exhaustive(&[&[]]);
    }

    #[test]
    fn single_slice() {
        assert_exhaustive(&[&[1, 2, 3, 4]]);
        assert_exhaustive(&[&[255, 128, 64]]);
    }

    #[test]
    fn two_slices() {
        assert_exhaustive(&[&[1, 2, 3], &[4, 5, 6]]);
        assert_exhaustive(&[&[0, 0, 0], &[0, 0, 1]]);
        assert_exhaustive(&[&[255, 255], &[255, 255]]);
    }

    #[test]
    fn three_slices() {
        assert_exhaustive(&[&[1, 2], &[3, 4], &[5, 6]]);
    }

    #[test]
    fn mixed_empty_slices() {
        assert_exhaustive(&[&[], &[1, 2, 3], &[]]);
        assert_exhaustive(&[&[1], &[], &[2]]);
    }

    #[test]
    fn unequal_lengths() {
        assert_exhaustive(&[&[1], &[2, 3, 4, 5]]);
        assert_exhaustive(&[&[1, 2, 3, 4], &[5]]);
    }

    #[test]
    fn matches_eval_coefficients_at() {
        let coefficients = [0x12u8, 0x34, 0x56, 0x78];

        for x in 0..=255u8 {
            assert_eq!(
                Polynomial::eval_coefficient_slices_at(&[&coefficients], x),
                Polynomial::eval_coefficients_at(&coefficients, x)
            );
        }
    }

    #[test]
    fn matches_eval_coefficient_iter_at() {
        let low = [1u8, 2, 3];
        let high = [4u8, 5, 6];

        for x in 0..=255u8 {
            assert_eq!(
                Polynomial::eval_coefficient_slices_at(&[&low, &high], x),
                Polynomial::eval_coefficient_iter_at(low.iter().chain(high.iter()).copied(), x)
            );
        }
    }
}
