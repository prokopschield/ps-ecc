use std::{array::IntoIter, iter::Take};

use crate::{polynomial::constants::POLYNOMIAL_MAX_COEFFICIENTS, Polynomial};

impl IntoIterator for Polynomial {
    type Item = u8;
    type IntoIter = Take<IntoIter<u8, POLYNOMIAL_MAX_COEFFICIENTS>>;

    fn into_iter(self) -> Self::IntoIter {
        self.coefficients
            .into_iter()
            .take(self.degree.saturating_add(1) as usize)
    }
}

impl<'a> IntoIterator for &'a Polynomial {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn iter_ref_polynomial() {
        let poly = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let coeffs: Vec<u8> = (&poly).into_iter().copied().collect();

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
    fn iter_with_adapters() {
        let poly = Polynomial::try_from(&[1u8, 2, 3, 4][..]).expect("valid polynomial");
        let doubled: Vec<u8> = (&poly).into_iter().map(|&c| c.wrapping_mul(2)).collect();

        assert_eq!(doubled, vec![2, 4, 6, 8]);
    }

    #[test]
    fn iter_empty_polynomial() {
        let poly = Polynomial::default();
        let coeffs: Vec<u8> = (&poly).into_iter().copied().collect();

        assert_eq!(coeffs, vec![0]);
    }

    #[test]
    fn iter_method() {
        let poly = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let coeffs: Vec<u8> = poly.iter().copied().collect();

        assert_eq!(coeffs, vec![1, 2, 3]);
    }

    #[test]
    fn into_iter_owned() {
        let poly = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let coeffs: Vec<u8> = poly.into_iter().collect();

        assert_eq!(coeffs, vec![1, 2, 3]);
    }

    #[test]
    fn into_iter_owned_in_for_loop() {
        let poly = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");
        let mut sum = 0u8;

        for coef in poly {
            sum ^= coef;
        }

        assert_eq!(sum, 4 ^ 5 ^ 6);
    }

    #[test]
    fn into_iter_exact_size() {
        let poly = Polynomial::try_from(&[1u8, 2, 3, 4, 5][..]).expect("valid polynomial");
        let iter = poly.into_iter();

        assert_eq!(iter.len(), 5);
    }

    #[test]
    fn into_iter_max_degree() {
        let coeffs = [1u8; 255];
        let poly = Polynomial::try_from(&coeffs[..]).expect("valid polynomial");

        assert_eq!(poly.into_iter().count(), 255);
    }

    #[test]
    fn into_iter_rev() {
        let poly = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let reversed: Vec<u8> = poly.into_iter().rev().collect();

        assert_eq!(reversed, vec![3, 2, 1]);
    }

    #[test]
    fn into_iter_matches_iter() {
        let poly = Polynomial::try_from(&[1u8, 2, 3, 4, 5][..]).expect("valid polynomial");
        let from_iter: Vec<u8> = poly.iter().copied().collect();
        let from_into_iter: Vec<u8> = poly.into_iter().collect();

        assert_eq!(from_iter, from_into_iter);
    }
}
