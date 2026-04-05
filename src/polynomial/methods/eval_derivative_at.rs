use crate::Polynomial;

impl Polynomial {
    /// Evaluates the formal derivative of this polynomial at the given point.
    ///
    /// In characteristic 2, only odd-degree terms contribute to the derivative.
    #[must_use]
    pub fn eval_derivative_at(&self, x: u8) -> u8 {
        Self::eval_coefficients_derivative_at(self.coefficients(), x)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        finite_field::{add, mul},
        Polynomial,
    };

    /// Asserts that `eval_derivative_at` matches `eval_coefficients_derivative_at` for all x.
    fn assert_exhaustive(p: &Polynomial) {
        for x in 0..=255u8 {
            assert_eq!(
                p.eval_derivative_at(x),
                Polynomial::eval_coefficients_derivative_at(p.coefficients(), x)
            );
        }
    }

    #[test]
    fn zero_polynomial() {
        assert_exhaustive(&Polynomial::default());
    }

    #[test]
    fn constant_polynomial() {
        let mut p = Polynomial::default();

        p.set(0, 42);

        assert_exhaustive(&p);
    }

    #[test]
    fn linear_polynomial() {
        let mut p = Polynomial::default();

        p.set(0, 5);
        p.set(1, 3);

        assert_exhaustive(&p);
    }

    #[test]
    fn quadratic_polynomial() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(1, 7);
        p.set(2, 1);

        assert_exhaustive(&p);
    }

    #[test]
    fn cubic_polynomial() {
        let mut p = Polynomial::default();

        p.set(0, 0x12);
        p.set(1, 0x34);
        p.set(2, 0x56);
        p.set(3, 0x78);

        assert_exhaustive(&p);
    }

    #[test]
    fn high_degree_polynomial() {
        let mut p = Polynomial::default();

        for i in 0..8 {
            p.set(i, i + 1);
        }

        assert_exhaustive(&p);
    }

    #[test]
    fn max_degree_polynomial() {
        let mut p = Polynomial::default();

        p.set(Polynomial::MAX_DEGREE, 1);

        assert_exhaustive(&p);
    }

    #[test]
    fn only_even_powers() {
        // p(x) = 1 + x² + x⁴, p'(x) = 0
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(2, 1);
        p.set(4, 1);

        for x in 0..=255u8 {
            assert_eq!(p.eval_derivative_at(x), 0);
        }
    }

    #[test]
    fn only_odd_powers() {
        // p(x) = x + x³, p'(x) = 1 + x²
        let mut p = Polynomial::default();

        p.set(1, 1);
        p.set(3, 1);

        for x in 0..=255u8 {
            let x2 = mul(x, x);
            let expected = add(1, x2);

            assert_eq!(p.eval_derivative_at(x), expected);
        }
    }

    #[test]
    fn derivative_computed_correctly() {
        // p(x) = 0x12 + 0x34*x + 0x56*x² + 0x78*x³ + 0x9A*x⁴
        // p'(x) = 0x34 + 0x78*x² (only odd-degree terms contribute)
        let mut p = Polynomial::default();

        p.set(0, 0x12);
        p.set(1, 0x34);
        p.set(2, 0x56);
        p.set(3, 0x78);
        p.set(4, 0x9A);

        for x in 0..=255u8 {
            let x2 = mul(x, x);
            let expected = add(0x34, mul(0x78, x2));

            assert_eq!(p.eval_derivative_at(x), expected);
        }
    }
}
