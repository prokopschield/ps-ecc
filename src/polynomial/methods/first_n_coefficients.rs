use crate::{polynomial::constants::POLYNOMIAL_MAX_COEFFICIENTS, Polynomial};

impl Polynomial {
    /// Returns the first `n` coefficients.
    ///
    /// The returned slice has length `min(n, 255)`.
    #[must_use]
    pub fn first_n_coefficients(&self, n: usize) -> &[u8] {
        &self.coefficients[..n.min(POLYNOMIAL_MAX_COEFFICIENTS)]
    }
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn exact_size() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(1, 2);
        p.set(2, 3);

        assert_eq!(p.first_n_coefficients(3), &[1, 2, 3]);
    }

    #[test]
    fn zeros_beyond_degree() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(1, 2);

        assert_eq!(p.first_n_coefficients(5), &[1, 2, 0, 0, 0]);
    }

    #[test]
    fn truncates() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(1, 2);
        p.set(2, 3);
        p.set(3, 4);
        p.set(4, 5);

        assert_eq!(p.first_n_coefficients(3), &[1, 2, 3]);
    }

    #[test]
    fn zero_polynomial() {
        let p = Polynomial::default();

        assert_eq!(p.first_n_coefficients(3), &[0, 0, 0]);
    }

    #[test]
    fn n_is_zero() {
        let mut p = Polynomial::default();

        p.set(0, 1);

        assert_eq!(p.first_n_coefficients(0), &[]);
    }

    #[test]
    fn n_is_255() {
        let mut p = Polynomial::default();

        p.set(0, 1);

        let result = p.first_n_coefficients(255);

        assert_eq!(result.len(), 255);
        assert_eq!(result[0], 1);
        assert!(result[1..].iter().all(|&x| x == 0));
    }

    #[test]
    fn n_exceeds_255_clamps() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(1, 2);

        let result = p.first_n_coefficients(300);

        assert_eq!(result.len(), 255);
    }
}
