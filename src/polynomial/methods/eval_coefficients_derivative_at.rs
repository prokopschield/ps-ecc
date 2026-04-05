use crate::{
    finite_field::{add, mul},
    Polynomial,
};

impl Polynomial {
    /// Evaluates the formal derivative of a polynomial at the given point.
    ///
    /// In characteristic 2, only odd-degree terms contribute to the derivative:
    /// p'(x) = a₁ + a₃x² + a₅x⁴ + ...
    ///
    /// Coefficients are ordered from degree 0 (constant term) to the leading term.
    #[must_use]
    pub fn eval_coefficients_derivative_at(coefficients: &[u8], x: u8) -> u8 {
        let mut result = 0u8;
        let x_sq = mul(x, x);
        let mut x_pow = 1u8;

        for k in 0..=coefficients.len() / 2 {
            let idx = 2 * k + 1;

            if idx < coefficients.len() {
                result = add(result, mul(coefficients[idx], x_pow));
                x_pow = mul(x_pow, x_sq);
            }
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

    /// Computes expected derivative value using naive evaluation.
    fn naive_eval_deriv(coefficients: &[u8], x: u8) -> u8 {
        let mut result = 0u8;
        let x_sq = mul(x, x);
        let mut x_pow = 1u8;

        for k in 0..=coefficients.len() / 2 {
            let idx = 2 * k + 1;

            if idx < coefficients.len() {
                result = add(result, mul(coefficients[idx], x_pow));
                x_pow = mul(x_pow, x_sq);
            }
        }

        result
    }

    /// Asserts that `eval_coefficients_derivative_at` matches naive evaluation for all x.
    fn assert_exhaustive(coefficients: &[u8]) {
        for x in 0..=255u8 {
            assert_eq!(
                Polynomial::eval_coefficients_derivative_at(coefficients, x),
                naive_eval_deriv(coefficients, x)
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
    fn constant_polynomial_has_zero_derivative() {
        for c in [1u8, 42, 128, 255] {
            for x in 0..=255u8 {
                assert_eq!(Polynomial::eval_coefficients_derivative_at(&[c], x), 0);
            }
        }
    }

    #[test]
    fn linear_polynomial() {
        // p(x) = a + bx, p'(x) = b
        for x in 0..=255u8 {
            assert_eq!(Polynomial::eval_coefficients_derivative_at(&[5, 3], x), 3);
            assert_eq!(Polynomial::eval_coefficients_derivative_at(&[0, 1], x), 1);
            assert_eq!(
                Polynomial::eval_coefficients_derivative_at(&[255, 255], x),
                255
            );
        }
    }

    #[test]
    fn quadratic_polynomial() {
        // p(x) = a + bx + cx², p'(x) = b (c*2x = 0 in char 2)
        assert_exhaustive(&[1, 0, 1]);
        assert_exhaustive(&[0, 5, 1]);
        assert_exhaustive(&[255, 128, 64]);
    }

    #[test]
    fn cubic_polynomial() {
        // p(x) = a + bx + cx² + dx³, p'(x) = b + dx²
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
    fn only_even_powers() {
        // p(x) = 1 + x² + x⁴, p'(x) = 0
        for x in 0..=255u8 {
            assert_eq!(
                Polynomial::eval_coefficients_derivative_at(&[1, 0, 1, 0, 1], x),
                0
            );
        }
    }

    #[test]
    fn only_odd_powers() {
        // p(x) = x + x³, p'(x) = 1 + x²
        let coefficients = [0u8, 1, 0, 1];

        for x in 0..=255u8 {
            let x2 = mul(x, x);
            let expected = add(1, x2);

            assert_eq!(
                Polynomial::eval_coefficients_derivative_at(&coefficients, x),
                expected
            );
        }
    }

    #[test]
    fn max_degree_polynomial() {
        let mut coefficients = [0u8; Polynomial::MAX_COEFFICIENTS as usize];

        coefficients[Polynomial::MAX_DEGREE as usize] = 1;

        assert_exhaustive(&coefficients);
    }

    #[test]
    fn matches_naive_eval_deriv() {
        let coefficients = [0x12u8, 0x34, 0x56, 0x78, 0x9A];

        for x in 0..=255u8 {
            assert_eq!(
                Polynomial::eval_coefficients_derivative_at(&coefficients, x),
                naive_eval_deriv(&coefficients, x)
            );
        }
    }
}
