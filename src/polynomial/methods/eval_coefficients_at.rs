use crate::{
    finite_field::{add, mul},
    Polynomial,
};

impl Polynomial {
    /// Evaluates a polynomial given its coefficients at the given point using Horner's method.
    ///
    /// Coefficients are ordered from degree 0 (constant term) to the leading term.
    #[must_use]
    pub fn eval_coefficients_at(coefficients: &[u8], x: u8) -> u8 {
        let mut result = 0u8;

        for &coef in coefficients.iter().rev() {
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
    fn naive_eval(coefficients: &[u8], x: u8) -> u8 {
        let mut result = 0u8;
        let mut x_pow = 1u8;

        for &c in coefficients {
            result = add(result, mul(c, x_pow));
            x_pow = mul(x_pow, x);
        }

        result
    }

    /// Asserts that `eval_coefficients_at` matches naive evaluation for all x in 0..=255.
    fn assert_exhaustive(coefficients: &[u8]) {
        for x in 0..=255u8 {
            assert_eq!(
                Polynomial::eval_coefficients_at(coefficients, x),
                naive_eval(coefficients, x)
            );
        }
    }

    #[test]
    fn empty_slice() {
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
}
