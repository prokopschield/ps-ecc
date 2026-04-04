use crate::Polynomial;

impl Polynomial {
    /// Evaluates this polynomial at the given point.
    #[must_use]
    pub fn eval_at(&self, x: u8) -> u8 {
        Self::eval_coefficients_at(self.coefficients(), x)
    }
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    /// Asserts that `eval_at` matches `eval_coefficients_at` for all x in 0..=255.
    fn assert_exhaustive(p: &Polynomial) {
        for x in 0..=255u8 {
            assert_eq!(
                p.eval_at(x),
                Polynomial::eval_coefficients_at(p.coefficients(), x)
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

        let mut p = Polynomial::default();

        p.set(0, 255);

        assert_exhaustive(&p);

        let mut p = Polynomial::default();

        p.set(0, 1);

        assert_exhaustive(&p);
    }

    #[test]
    fn linear_polynomial() {
        let mut p = Polynomial::default();

        p.set(0, 5);
        p.set(1, 3);

        assert_exhaustive(&p);

        let mut p = Polynomial::default();

        p.set(0, 0);
        p.set(1, 1);

        assert_exhaustive(&p);

        let mut p = Polynomial::default();

        p.set(0, 255);
        p.set(1, 255);

        assert_exhaustive(&p);
    }

    #[test]
    fn quadratic_polynomial() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(1, 0);
        p.set(2, 1);

        assert_exhaustive(&p);

        let mut p = Polynomial::default();

        p.set(0, 0);
        p.set(1, 0);
        p.set(2, 1);

        assert_exhaustive(&p);

        let mut p = Polynomial::default();

        p.set(0, 255);
        p.set(1, 128);
        p.set(2, 64);

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

        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(1, 1);
        p.set(2, 1);
        p.set(3, 1);

        assert_exhaustive(&p);

        let mut p = Polynomial::default();

        p.set(3, 1);

        assert_exhaustive(&p);
    }

    #[test]
    fn high_degree_polynomial() {
        let mut p = Polynomial::default();

        for i in 0..8 {
            p.set(i, i + 1);
        }

        assert_exhaustive(&p);

        let mut p = Polynomial::default();

        p.set(7, 1);

        assert_exhaustive(&p);
    }

    #[test]
    fn sparse_polynomial() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(4, 1);

        assert_exhaustive(&p);

        let mut p = Polynomial::default();

        p.set(1, 1);
        p.set(3, 1);

        assert_exhaustive(&p);
    }

    #[test]
    fn all_coefficients_max() {
        let mut p = Polynomial::default();

        p.set(0, 255);
        p.set(1, 255);
        p.set(2, 255);
        p.set(3, 255);

        assert_exhaustive(&p);
    }

    #[test]
    fn max_degree_monomial() {
        let mut p = Polynomial::default();

        p.set(Polynomial::MAX_DEGREE, 1);

        assert_exhaustive(&p);
    }

    #[test]
    fn max_degree_all_ones() {
        let mut p = Polynomial::default();

        for i in 0..=Polynomial::MAX_DEGREE {
            p.set(i, 1);
        }

        assert_exhaustive(&p);
    }

    #[test]
    fn max_degree_all_max() {
        let mut p = Polynomial::default();

        for i in 0..=Polynomial::MAX_DEGREE {
            p.set(i, 255);
        }

        assert_exhaustive(&p);
    }
}
