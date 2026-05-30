use crate::Polynomial;

impl Polynomial {
    /// Returns the coefficients of this polynomial, from degree 0 to the leading term.
    ///
    /// # Panics
    ///
    /// Panics on polynomials of degree 255, which are invalid in GF(256).
    #[must_use]
    pub fn coefficients(&self) -> &[u8] {
        #[allow(clippy::expect_used)]
        self.coefficients.get(..self.degree() as usize + 1).expect(
            "BUG: Polynomial::coefficients called on an invalid GF(256) polynomial of degree 255",
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn default_polynomial_has_one_coefficient() {
        let p = Polynomial::default();

        assert_eq!(p.coefficients(), &[0]);
    }

    #[test]
    fn coefficients_match_set_values() {
        let mut p = Polynomial::default();

        p.set(0, 5);
        p.set(1, 3);
        p.set(2, 7);

        assert_eq!(p.coefficients(), &[5, 3, 7]);
    }

    #[test]
    fn sparse_polynomial_includes_zeros() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(3, 2);

        assert_eq!(p.coefficients(), &[1, 0, 0, 2]);
    }

    #[test]
    fn coefficients_length_is_degree_plus_one() {
        let mut p = Polynomial::default();

        p.set(5, 1);

        assert_eq!(p.coefficients().len(), 6);
    }

    #[test]
    fn coefficients_after_trim() {
        let mut p = Polynomial::default();

        p.set(4, 1);
        p.set(4, 0);

        assert_eq!(p.coefficients(), &[0]);
    }
}
