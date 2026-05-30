use crate::Polynomial;

impl Polynomial {
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.degree
    }
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn zero_polynomial() {
        let p = Polynomial::default();

        assert_eq!(p.degree(), 0);
    }

    #[test]
    fn constant_polynomial() {
        let mut p = Polynomial::default();

        p.set(0, 42);

        assert_eq!(p.degree(), 0);
    }

    #[test]
    fn linear() {
        let mut p = Polynomial::default();

        p.set(1, 1);

        assert_eq!(p.degree(), 1);
    }

    #[test]
    fn quadratic() {
        let mut p = Polynomial::default();

        p.set(2, 1);

        assert_eq!(p.degree(), 2);
    }

    #[test]
    fn max_degree() {
        let mut p = Polynomial::default();

        p.set(254, 1);

        assert_eq!(p.degree(), 254);
    }

    #[test]
    fn degree_after_trim() {
        let mut p = Polynomial::default();

        p.set(4, 1);
        p.set(4, 0);

        assert_eq!(p.degree(), 0);
    }

    #[test]
    fn degree_after_partial_trim() {
        let mut p = Polynomial::default();

        p.set(2, 5);
        p.set(5, 3);
        p.set(5, 0);

        assert_eq!(p.degree(), 2);
    }

    #[test]
    fn degree_increases() {
        let mut p = Polynomial::default();

        p.set(1, 1);

        assert_eq!(p.degree(), 1);

        p.set(3, 1);

        assert_eq!(p.degree(), 3);

        p.set(7, 1);

        assert_eq!(p.degree(), 7);
    }

    #[test]
    fn degree_unaffected_by_lower_sets() {
        let mut p = Polynomial::default();

        p.set(5, 1);
        p.set(2, 9);
        p.set(0, 3);

        assert_eq!(p.degree(), 5);
    }

    #[test]
    fn degree_with_sparse_coefficients() {
        let mut p = Polynomial::default();

        p.set(0, 1);
        p.set(10, 1);
        p.set(5, 0);

        assert_eq!(p.degree(), 10);
    }
}
