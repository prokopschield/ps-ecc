use crate::{
    finite_field::{add, mul},
    Polynomial,
};

impl Polynomial {
    /// Evaluates a polynomial from an iterator of coefficients at the given point.
    ///
    /// Coefficients are yielded from degree 0 (constant term) to the leading term.
    /// Requires `DoubleEndedIterator` for Horner's method.
    #[must_use]
    pub fn eval_coefficient_iter_at<I>(coefficients: I, x: u8) -> u8
    where
        I: DoubleEndedIterator<Item = u8>,
    {
        let mut result = 0u8;

        for coef in coefficients.rev() {
            result = add(mul(result, x), coef);
        }

        result
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

    /// Asserts that `eval_coefficient_iter_at` matches naive evaluation for all x in 0..=255.
    fn assert_exhaustive(coefficients: &[u8]) {
        for x in 0..=255u8 {
            assert_eq!(
                Polynomial::eval_coefficient_iter_at(coefficients.iter().copied(), x),
                naive_eval(coefficients.iter().copied(), x)
            );
        }
    }

    /// Asserts that `eval_coefficient_iter_at` with chained slices matches naive evaluation for all x.
    fn assert_exhaustive_chained(low: &[u8], high: &[u8]) {
        for x in 0..=255u8 {
            let iter = low.iter().chain(high.iter()).copied();

            assert_eq!(
                Polynomial::eval_coefficient_iter_at(low.iter().chain(high.iter()).copied(), x),
                naive_eval(iter, x)
            );
        }
    }

    #[test]
    fn empty_iterator() {
        assert_exhaustive(&[]);
    }

    #[test]
    fn zero_polynomial() {
        assert_exhaustive(&[0]);
    }

    #[test]
    fn constant_polynomial() {
        assert_exhaustive(&[1]);
        assert_exhaustive(&[42]);
        assert_exhaustive(&[128]);
        assert_exhaustive(&[255]);
    }

    #[test]
    fn linear_polynomial() {
        assert_exhaustive(&[5, 3]);
        assert_exhaustive(&[0, 1]);
        assert_exhaustive(&[255, 255]);
    }

    #[test]
    fn quadratic_polynomial() {
        assert_exhaustive(&[1, 0, 1]);
        assert_exhaustive(&[0, 0, 1]);
        assert_exhaustive(&[255, 128, 64]);
    }

    #[test]
    fn cubic_polynomial() {
        assert_exhaustive(&[0x12, 0x34, 0x56, 0x78]);
        assert_exhaustive(&[1, 1, 1, 1]);
        assert_exhaustive(&[0, 0, 0, 1]);
    }

    #[test]
    fn high_degree_polynomial() {
        assert_exhaustive(&[1, 2, 3, 4, 5, 6, 7, 8]);
        assert_exhaustive(&[0, 0, 0, 0, 0, 0, 0, 1]);
    }

    #[test]
    fn sparse_polynomial() {
        assert_exhaustive(&[1, 0, 0, 0, 1]);
        assert_exhaustive(&[0, 1, 0, 1, 0]);
    }

    #[test]
    fn all_coefficients_max() {
        assert_exhaustive(&[255, 255, 255, 255]);
    }

    #[test]
    fn max_degree_monomial() {
        let mut coefficients = [0u8; Polynomial::MAX_COEFFICIENTS as usize];
        coefficients[Polynomial::MAX_DEGREE as usize] = 1;

        assert_exhaustive(&coefficients);
    }

    #[test]
    fn max_degree_all_ones() {
        let coefficients = [1u8; Polynomial::MAX_COEFFICIENTS as usize];

        assert_exhaustive(&coefficients);
    }

    #[test]
    fn max_degree_all_max() {
        let coefficients = [255u8; Polynomial::MAX_COEFFICIENTS as usize];

        assert_exhaustive(&coefficients);
    }

    #[test]
    fn chained_empty_slices() {
        assert_exhaustive_chained(&[], &[]);
    }

    #[test]
    fn chained_one_empty() {
        assert_exhaustive_chained(&[1, 2, 3], &[]);
        assert_exhaustive_chained(&[], &[1, 2, 3]);
    }

    #[test]
    fn chained_equal_length() {
        assert_exhaustive_chained(&[1, 2, 3], &[4, 5, 6]);
        assert_exhaustive_chained(&[0, 0, 0], &[0, 0, 1]);
        assert_exhaustive_chained(&[255, 255], &[255, 255]);
    }

    #[test]
    fn chained_unequal_length() {
        assert_exhaustive_chained(&[1], &[2, 3, 4, 5]);
        assert_exhaustive_chained(&[1, 2, 3, 4], &[5]);
    }

    #[test]
    fn matches_eval_coefficients_at() {
        let coefficients = [0x12u8, 0x34, 0x56, 0x78];

        for x in 0..=255u8 {
            assert_eq!(
                Polynomial::eval_coefficient_iter_at(coefficients.iter().copied(), x),
                Polynomial::eval_coefficients_at(&coefficients, x)
            );
        }
    }
}
