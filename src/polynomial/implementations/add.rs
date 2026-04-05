use std::{
    borrow::Borrow,
    ops::{Add, BitXor},
};

use crate::{error::PolynomialXorError, Polynomial};

impl<B: Borrow<u8>, I: IntoIterator<Item = B>> Add<I> for Polynomial {
    type Output = Result<Self, PolynomialXorError>;

    /// Adds coefficients from an iterator to this polynomial.
    ///
    /// In GF(256), addition is XOR. This delegates to `BitXor`.
    ///
    /// # Errors
    ///
    /// Returns `TooManyCoefficients` if the iterator yields more than 255 elements.
    fn add(self, rhs: I) -> Self::Output {
        BitXor::bitxor(self, rhs)
    }
}

impl<B: Borrow<u8>, I: IntoIterator<Item = B>> Add<I> for &Polynomial {
    type Output = Result<Polynomial, PolynomialXorError>;

    fn add(self, rhs: I) -> Self::Output {
        BitXor::bitxor(*self, rhs)
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::ops::Add;

    use crate::{error::PolynomialXorError, Polynomial};

    #[test]
    fn add_two_polynomials() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        let result = a.add(b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn add_different_degrees() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[3u8, 4, 5, 6][..]).expect("valid polynomial");

        let result = a.add(b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 3, 2 ^ 4, 5, 6]);
    }

    #[test]
    fn add_cancels_to_zero() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let result = a.add(b).expect("polynomials are bounded");

        assert_eq!(result.degree(), 0);
        assert_eq!(result.coefficients(), &[0]);
    }

    #[test]
    fn add_with_references() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        let result = (&a).add(&b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn add_slice_via_iter() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let slice: &[u8] = &[4, 5, 6];

        let result = a.add(slice.iter()).expect("valid slice");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn add_slice_different_lengths() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let slice: &[u8] = &[3, 4, 5, 6];

        let result = a.add(slice.iter()).expect("valid slice");

        assert_eq!(result.coefficients(), &[1 ^ 3, 2 ^ 4, 5, 6]);
    }

    #[test]
    fn add_iterator() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let iter = [4u8, 5, 6].into_iter();

        let result = a.add(iter).expect("valid iterator");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn add_vec() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let vec = vec![4u8, 5, 6, 7];

        let result = a.add(vec).expect("valid vec");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6, 7]);
    }

    #[test]
    fn add_iterator_too_long_returns_error() {
        let a = Polynomial::default();
        let iter = std::iter::repeat_n(1u8, 256);

        let result = a.add(iter);

        assert_eq!(result, Err(PolynomialXorError::TooManyCoefficients));
    }
}
