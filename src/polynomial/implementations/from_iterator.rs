use crate::Polynomial;

#[allow(clippy::cast_possible_truncation)] // len(coefficients) < u8::MAX
impl FromIterator<u8> for Polynomial {
    /// Collects coefficients from degree 0 (constant term) to the leading term.
    ///
    /// Consumes at most 255 elements; excess elements are ignored.
    /// Leading zeros are trimmed. Safe for infinite iterators.
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut poly = Self::default();

        poly.degree = poly
            .coefficients
            .iter_mut()
            .zip(iter)
            .map(|(coef, value)| *coef = value)
            .count() as u8;

        poly.trim_degree();

        poly
    }
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn empty_iterator() {
        let poly: Polynomial = std::iter::empty().collect();

        assert_eq!(poly.degree(), 0);
        assert_eq!(poly.coefficients(), &[0]);
    }

    #[test]
    fn collect_vec() {
        let poly: Polynomial = vec![1u8, 2, 3].into_iter().collect();

        assert_eq!(poly.degree(), 2);
        assert_eq!(poly.coefficients(), &[1, 2, 3]);
    }

    #[test]
    fn collect_slice_copied() {
        let slice: &[u8] = &[1, 2, 3];
        let poly: Polynomial = slice.iter().copied().collect();

        assert_eq!(poly.coefficients(), &[1, 2, 3]);
    }

    #[test]
    fn trims_leading_zeros() {
        let poly: Polynomial = [1u8, 2, 0, 0].into_iter().collect();

        assert_eq!(poly.degree(), 1);
        assert_eq!(poly.coefficients(), &[1, 2]);
    }

    #[test]
    fn max_size() {
        let poly: Polynomial = [1u8; 255].into_iter().collect();

        assert_eq!(poly.degree(), 254);
        assert_eq!(poly.coefficients().len(), 255);
    }

    #[test]
    fn truncates_excess_elements() {
        let poly: Polynomial = [1u8; 256].into_iter().collect();

        assert_eq!(poly.degree(), 254);
        assert_eq!(poly.coefficients().len(), 255);
    }

    #[test]
    fn infinite_iterator() {
        let poly: Polynomial = std::iter::repeat(1u8).collect();

        assert_eq!(poly.degree(), 254);
        assert_eq!(poly.coefficients().len(), 255);
    }
}
