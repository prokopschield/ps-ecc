use std::ops::Div;

use crate::{error::PolynomialDivError, Polynomial};

impl<D: AsRef<[u8]>> Div<D> for &Polynomial {
    type Output = Result<Polynomial, PolynomialDivError>;

    /// Divides two polynomials in GF(256), returning the quotient.
    ///
    /// # Errors
    ///
    /// Returns `ZeroDivisor` if the divisor is the zero polynomial.
    fn div(self, rhs: D) -> Self::Output {
        let (quot, _rem) = self.div_rem(rhs)?;

        Ok(quot)
    }
}

impl<D: AsRef<[u8]>> Div<D> for Polynomial {
    type Output = Result<Self, PolynomialDivError>;

    fn div(self, rhs: D) -> Self::Output {
        let (quot, _rem) = self.div_rem(rhs)?;

        Ok(quot)
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::ops::{Div, Mul};

    use crate::{error::PolynomialDivError, Polynomial};

    #[test]
    fn div_returns_quotient() {
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        let quot = a.div(&b).expect("division succeeds");

        assert_eq!(quot.coefficients(), &[1, 1]);
    }

    #[test]
    fn div_by_slice() {
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");
        let b: &[u8] = &[1, 1];

        let quot = a.div(b).expect("division succeeds");

        assert_eq!(quot.coefficients(), &[1, 1]);
    }

    #[test]
    fn div_by_array() {
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");

        let quot = a.div([1u8, 1]).expect("division succeeds");

        assert_eq!(quot.coefficients(), &[1, 1]);
    }

    #[test]
    fn div_by_zero() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let zero = Polynomial::default();

        let result = a.div(&zero);

        assert_eq!(result, Err(PolynomialDivError::ZeroDivisor));
    }

    #[test]
    fn div_mul_roundtrip() {
        let a = Polynomial::try_from(&[3u8, 5, 7][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[2u8, 4][..]).expect("valid polynomial");

        let product = (&a).mul(&b).expect("multiplication succeeds");
        let quot = product.div(&b).expect("division succeeds");

        assert_eq!(quot, a);
    }

    #[test]
    fn div_ref_lhs() {
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        let quot = (&a).div(&b).expect("division succeeds");

        assert_eq!(quot.coefficients(), &[1, 1]);
        // a is still usable
        assert_eq!(a.coefficients(), &[1, 0, 1]);
    }
}
