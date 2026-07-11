use std::ops::Rem;

use crate::{error::PolynomialDivError, Polynomial};

impl<D: AsRef<[u8]>> Rem<D> for &Polynomial {
    type Output = Result<Polynomial, PolynomialDivError>;

    /// Computes the remainder of polynomial division over GF(256).
    ///
    /// # Errors
    ///
    /// Returns `ZeroDivisor` if the divisor is the zero polynomial.
    fn rem(self, rhs: D) -> Self::Output {
        let (_quot, rem) = self.div_rem(rhs)?;

        Ok(rem)
    }
}

impl<D: AsRef<[u8]>> Rem<D> for Polynomial {
    type Output = Result<Self, PolynomialDivError>;

    fn rem(self, rhs: D) -> Self::Output {
        let (_quot, rem) = self.div_rem(rhs)?;

        Ok(rem)
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::ops::{Div, Mul, Rem};

    use crate::{error::PolynomialDivError, Polynomial};

    #[test]
    fn rem_returns_remainder() {
        // (x^2 + x + 1) % (x + 1) has a non-zero remainder
        let a = Polynomial::try_from(&[1u8, 1, 1][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        let rem = a.rem(&b).expect("division succeeds");

        // Verify via div_rem
        let (_quot, expected_rem) = a.div_rem(b).expect("division succeeds");

        assert_eq!(rem, expected_rem);
    }

    #[test]
    fn rem_zero_when_divisible() {
        // (x^2 + 1) % (x + 1) = 0 in GF(256)
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        let rem = a.rem(&b).expect("division succeeds");

        assert_eq!(rem.coefficients(), &[0]);
    }

    #[test]
    fn rem_by_slice() {
        let a = Polynomial::try_from(&[1u8, 1, 1][..]).expect("valid polynomial");
        let b: &[u8] = &[1, 1];

        let rem = a.rem(b).expect("division succeeds");
        let (_quot, expected_rem) = a.div_rem(b).expect("division succeeds");

        assert_eq!(rem, expected_rem);
    }

    #[test]
    fn rem_by_array() {
        let a = Polynomial::try_from(&[1u8, 1, 1][..]).expect("valid polynomial");

        let rem = a.rem([1u8, 1]).expect("division succeeds");
        let (_quot, expected_rem) = a.div_rem([1u8, 1]).expect("division succeeds");

        assert_eq!(rem, expected_rem);
    }

    #[test]
    fn rem_by_zero() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let zero = Polynomial::default();

        let result = a.rem(&zero);

        assert_eq!(result, Err(PolynomialDivError::ZeroDivisor));
    }

    #[test]
    fn rem_ref_lhs() {
        let a = Polynomial::try_from(&[1u8, 1, 1][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        let rem = (&a).rem(&b).expect("division succeeds");
        let (_quot, expected_rem) = a.div_rem(b).expect("division succeeds");

        assert_eq!(rem, expected_rem);
    }

    #[test]
    fn rem_smaller_than_divisor() {
        // When dividend degree < divisor degree, remainder = dividend
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let rem = a.rem(&b).expect("division succeeds");

        assert_eq!(rem.coefficients(), &[1, 2]);
    }

    #[test]
    fn div_rem_identity() {
        // Verify: (a / b) * b + (a % b) = a
        let a = Polynomial::try_from(&[3u8, 5, 7, 11][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[2u8, 4][..]).expect("valid polynomial");

        let quot = (&a).div(&b).expect("division succeeds");
        let rem = (&a).rem(&b).expect("division succeeds");

        let product = quot.mul(&b).expect("multiplication succeeds");
        let reconstructed = (product + rem).expect("addition succeeds");

        assert_eq!(reconstructed, a);
    }
}
