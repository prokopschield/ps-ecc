use std::{
    borrow::Borrow,
    ops::{BitXor, Sub},
};

use crate::{error::PolynomialXorError, Polynomial};

impl<B: Borrow<u8>, I: IntoIterator<Item = B>> Sub<I> for Polynomial {
    type Output = Result<Self, PolynomialXorError>;

    /// Subtracts coefficients from an iterator from this polynomial.
    ///
    /// In GF(256), subtraction is equivalent to addition, which is XOR.
    /// This delegates to `BitXor`.
    ///
    /// # Errors
    ///
    /// Returns `TooManyCoefficients` if the iterator yields more than 255 elements.
    fn sub(self, rhs: I) -> Self::Output {
        BitXor::bitxor(self, rhs)
    }
}

impl<B: Borrow<u8>, I: IntoIterator<Item = B>> Sub<I> for &Polynomial {
    type Output = Result<Polynomial, PolynomialXorError>;

    fn sub(self, rhs: I) -> Self::Output {
        BitXor::bitxor(*self, rhs)
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::ops::Sub;

    use crate::{error::PolynomialXorError, Polynomial};

    #[test]
    fn sub_two_polynomials() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        let result = a.sub(b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn sub_different_degrees() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[3u8, 4, 5, 6][..]).expect("valid polynomial");

        let result = a.sub(b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 3, 2 ^ 4, 5, 6]);
    }

    #[test]
    fn sub_cancels_to_zero() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let result = a.sub(b).expect("polynomials are bounded");

        assert_eq!(result.degree(), 0);
        assert_eq!(result.coefficients(), &[0]);
    }

    #[test]
    fn sub_with_references() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        let result = (&a).sub(&b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn sub_slice_via_iter() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let slice: &[u8] = &[4, 5, 6];

        let result = a.sub(slice.iter()).expect("valid slice");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn sub_slice_different_lengths() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let slice: &[u8] = &[3, 4, 5, 6];

        let result = a.sub(slice.iter()).expect("valid slice");

        assert_eq!(result.coefficients(), &[1 ^ 3, 2 ^ 4, 5, 6]);
    }

    #[test]
    fn sub_iterator() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let iter = [4u8, 5, 6].into_iter();

        let result = a.sub(iter).expect("valid iterator");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn sub_vec() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let vec = vec![4u8, 5, 6, 7];

        let result = a.sub(vec).expect("valid vec");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6, 7]);
    }

    #[test]
    fn sub_iterator_too_long_returns_error() {
        let a = Polynomial::default();
        let iter = std::iter::repeat_n(1u8, 256);

        let result = a.sub(iter);

        assert_eq!(result, Err(PolynomialXorError::TooManyCoefficients));
    }

    #[test]
    fn sub_equals_add() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        let sub_result = a.sub(&b).expect("polynomials are bounded");

        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let add_result = std::ops::Add::add(a, &b).expect("polynomials are bounded");

        assert_eq!(sub_result, add_result);
    }
}
