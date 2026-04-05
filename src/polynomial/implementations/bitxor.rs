use std::{borrow::Borrow, ops::BitXor};

use crate::{error::PolynomialXorError, finite_field::add, Polynomial};

impl<B: Borrow<u8>, I: IntoIterator<Item = B>> BitXor<I> for Polynomial {
    type Output = Result<Self, PolynomialXorError>;

    /// XORs coefficients from an iterator into this polynomial.
    ///
    /// # Errors
    ///
    /// Returns `TooManyCoefficients` if the iterator yields more than 255 elements.
    fn bitxor(self, rhs: I) -> Self::Output {
        let mut result = self;

        for (idx, coef) in rhs.into_iter().enumerate() {
            if idx >= Self::MAX_COEFFICIENTS as usize {
                return Err(PolynomialXorError::TooManyCoefficients);
            }

            result.coefficients[idx] = add(result.coefficients[idx], *coef.borrow());

            #[allow(clippy::cast_possible_truncation)]
            let idx_u8 = idx as u8;

            result.degree = result.degree.max(idx_u8);
        }

        result.trim_degree();

        Ok(result)
    }
}

impl<B: Borrow<u8>, I: IntoIterator<Item = B>> BitXor<I> for &Polynomial {
    type Output = Result<Polynomial, PolynomialXorError>;

    fn bitxor(self, rhs: I) -> Self::Output {
        *self ^ rhs
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::ops::BitXor;

    use crate::{error::PolynomialXorError, Polynomial};

    #[test]
    fn xor_two_polynomials() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        let result = a.bitxor(b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn xor_different_degrees() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[3u8, 4, 5, 6][..]).expect("valid polynomial");

        let result = a.bitxor(b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 3, 2 ^ 4, 5, 6]);
    }

    #[test]
    fn xor_cancels_to_zero() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let result = a.bitxor(b).expect("polynomials are bounded");

        assert_eq!(result.degree(), 0);
        assert_eq!(result.coefficients(), &[0]);
    }

    #[test]
    fn xor_with_references() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        let result = (&a).bitxor(&b).expect("polynomials are bounded");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn xor_slice_via_iter() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let slice: &[u8] = &[4, 5, 6];

        let result = a.bitxor(slice.iter()).expect("valid slice");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn xor_vec() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let vec = vec![4u8, 5, 6, 7];

        let result = a.bitxor(vec).expect("valid vec");

        assert_eq!(result.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6, 7]);
    }

    #[test]
    fn xor_iterator_too_long_returns_error() {
        let a = Polynomial::default();
        let iter = std::iter::repeat_n(1u8, 256);

        let result = a.bitxor(iter);

        assert_eq!(result, Err(PolynomialXorError::TooManyCoefficients));
    }
}
