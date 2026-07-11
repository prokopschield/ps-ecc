use crate::Polynomial;

impl Polynomial {
    /// Returns an iterator over the coefficients, from degree 0 to the
    /// leading term.
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.coefficients().iter()
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn iter_basic() {
        let poly = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let coeffs: Vec<u8> = poly.iter().copied().collect();

        assert_eq!(coeffs, vec![1, 2, 3]);
    }

    #[test]
    fn iter_in_for_loop() {
        let poly = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");
        let mut sum = 0u8;

        for &coef in &poly {
            sum ^= coef;
        }

        assert_eq!(sum, 4 ^ 5 ^ 6);
    }

    #[test]
    fn iter_empty_polynomial() {
        let poly = Polynomial::default();
        let coeffs: Vec<u8> = poly.iter().copied().collect();

        assert_eq!(coeffs, vec![0]);
    }

    #[test]
    fn iter_with_adapters() {
        let poly = Polynomial::try_from(&[1u8, 2, 3, 4][..]).expect("valid polynomial");
        let doubled: Vec<u8> = poly.iter().map(|&c| c.wrapping_mul(2)).collect();

        assert_eq!(doubled, vec![2, 4, 6, 8]);
    }

    #[test]
    fn iter_length_matches_degree_plus_one() {
        let poly = Polynomial::try_from(&[1u8, 2, 3, 4, 5][..]).expect("valid polynomial");

        assert_eq!(poly.iter().count(), 5);
        assert_eq!(poly.iter().count(), poly.degree() as usize + 1);
    }

    #[test]
    fn iter_max_degree() {
        let coeffs = [1u8; 255];
        let poly = Polynomial::try_from(&coeffs[..]).expect("valid polynomial");

        assert_eq!(poly.iter().count(), 255);
    }

    #[test]
    fn iter_sparse_polynomial() {
        let poly = Polynomial::try_from(&[0u8, 0, 0, 1][..]).expect("valid polynomial");
        let coeffs: Vec<u8> = poly.iter().copied().collect();

        assert_eq!(coeffs, vec![0, 0, 0, 1]);
    }

    #[test]
    fn iter_is_borrowing() {
        let poly = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let _ = poly.iter();
        let _ = poly.iter();

        assert_eq!(poly.degree(), 2);
    }
}
