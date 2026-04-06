use crate::{
    error::PolynomialDivError,
    finite_field::{inv, mul, sub},
    Polynomial,
};

impl Polynomial {
    /// Divides this polynomial by a divisor, returning `(quotient, remainder)`.
    ///
    /// The divisor is interpreted as coefficients in ascending degree order.
    /// Trailing zeros are ignored when determining the divisor's degree.
    ///
    /// # Errors
    ///
    /// Returns `ZeroDivisor` if the divisor is the zero polynomial.
    pub fn div_rem(&self, divisor: impl AsRef<[u8]>) -> Result<(Self, Self), PolynomialDivError> {
        let divisor = divisor.as_ref();

        let Some(divisor_degree) = divisor.iter().rposition(|&c| c != 0) else {
            return Err(PolynomialDivError::ZeroDivisor);
        };

        let divisor_lc_inv = inv(divisor[divisor_degree])?.get();

        let Some(dividend_degree) = self.coefficients().iter().rposition(|&c| c != 0) else {
            return Ok((Self::default(), Self::default()));
        };

        if dividend_degree < divisor_degree {
            return Ok((Self::default(), *self));
        }

        let quot_degree = dividend_degree - divisor_degree;

        #[allow(clippy::cast_possible_truncation)]
        let mut quot = Self {
            degree: quot_degree as u8,
            ..Self::default()
        };

        let mut rem = *self;
        let mut deg = dividend_degree;

        while deg >= divisor_degree {
            let lead_coef = rem.coefficients[deg];

            if lead_coef == 0 {
                if deg == 0 {
                    break;
                }
                deg -= 1;
                continue;
            }

            let ratio = mul(lead_coef, divisor_lc_inv);
            let shift = deg - divisor_degree;

            quot.coefficients[shift] = ratio;
            rem.coefficients[deg] = 0;

            for (i, &divisor_coef) in divisor[..divisor_degree].iter().enumerate() {
                rem.coefficients[shift + i] =
                    sub(rem.coefficients[shift + i], mul(ratio, divisor_coef));
            }

            if deg == 0 {
                break;
            }
            deg -= 1;
        }

        quot.trim_degree();
        rem.trim_degree();

        Ok((quot, rem))
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::ops::Mul;

    use crate::{error::PolynomialDivError, Polynomial};

    #[test]
    fn div_rem_basic() {
        // (x^2 + 1) / (x + 1)
        // In GF(256), x + 1 divides x^2 + 1 with quotient x + 1 and remainder 0
        // because (x + 1)^2 = x^2 + 2x + 1 = x^2 + 1 (since 2 = 0 in GF(2))
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");

        let (quot, rem) = a.div_rem([1, 1]).expect("division succeeds");

        assert_eq!(quot.coefficients(), &[1, 1]);
        assert_eq!(rem.coefficients(), &[0]);
    }

    #[test]
    fn div_rem_by_one() {
        let a = Polynomial::try_from(&[5u8, 10, 15][..]).expect("valid polynomial");

        let (quot, rem) = a.div_rem([1]).expect("division succeeds");

        assert_eq!(quot.coefficients(), &[5, 10, 15]);
        assert_eq!(rem.coefficients(), &[0]);
    }

    #[test]
    fn div_rem_by_zero() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let result = a.div_rem([0]);

        assert_eq!(result, Err(PolynomialDivError::ZeroDivisor));
    }

    #[test]
    fn div_rem_by_empty_slice() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let result = a.div_rem([]);

        assert_eq!(result, Err(PolynomialDivError::ZeroDivisor));
    }

    #[test]
    fn div_rem_zero_by_nonzero() {
        let zero = Polynomial::default();

        let (quot, rem) = zero.div_rem([1, 2]).expect("division succeeds");

        assert_eq!(quot.degree(), 0);
        assert_eq!(quot.coefficients(), &[0]);
        assert_eq!(rem.coefficients(), &[0]);
    }

    #[test]
    fn div_rem_smaller_by_larger() {
        // Dividend degree < divisor degree => quotient is 0, remainder is dividend
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");

        let (quot, rem) = a.div_rem([1, 2, 3]).expect("division succeeds");

        assert_eq!(quot.degree(), 0);
        assert_eq!(quot.coefficients(), &[0]);
        assert_eq!(rem.coefficients(), &[1, 2]);
    }

    #[test]
    fn div_rem_mul_roundtrip() {
        // (a * b) / b = (a, 0) when b is not zero
        let a = Polynomial::try_from(&[3u8, 5, 7][..]).expect("valid polynomial");
        let b: &[u8] = &[2, 4];

        let product = (&a).mul(b).expect("multiplication succeeds");
        let (quot, rem) = product.div_rem(b).expect("division succeeds");

        assert_eq!(quot, a);
        assert_eq!(rem.coefficients(), &[0]);
    }

    #[test]
    fn div_rem_with_remainder() {
        // (x^2 + x + 1) / (x + 1) should have a remainder
        let a = Polynomial::try_from(&[1u8, 1, 1][..]).expect("valid polynomial");
        let b: &[u8] = &[1, 1];

        let (quot, rem) = a.div_rem(b).expect("division succeeds");

        // Verify: quot * b + rem = a
        let product = quot.mul(b).expect("multiplication succeeds");
        let reconstructed = (product + rem).expect("addition succeeds");

        assert_eq!(reconstructed, a);
    }

    #[test]
    fn div_rem_by_constant() {
        // Dividing by a constant should scale all coefficients, remainder 0
        let a = Polynomial::try_from(&[6u8, 12, 18][..]).expect("valid polynomial");

        let (quot, rem) = a.div_rem([2]).expect("division succeeds");

        assert_eq!(rem.coefficients(), &[0]);

        // In GF(256), division by 2 is multiplication by inverse of 2
        // We verify by checking that quot * 2 = a
        let check = quot.mul([2u8]).expect("multiplication succeeds");

        assert_eq!(check.coefficients(), &[6, 12, 18]);
    }

    #[test]
    fn div_rem_with_polynomial() {
        // Verify we can pass &Polynomial directly
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        let (quot, rem) = a.div_rem(b).expect("division succeeds");

        assert_eq!(quot.coefficients(), &[1, 1]);
        assert_eq!(rem.coefficients(), &[0]);
    }

    #[test]
    fn div_rem_trailing_zeros_ignored() {
        // Trailing zeros in divisor should be ignored
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");

        let (quot1, rem1) = a.div_rem([1, 1]).expect("division succeeds");
        let (quot2, rem2) = a.div_rem([1, 1, 0, 0, 0]).expect("division succeeds");

        assert_eq!(quot1, quot2);
        assert_eq!(rem1, rem2);
    }

    #[test]
    fn div_rem_same_degree() {
        // Dividend and divisor have the same degree
        let a = Polynomial::try_from(&[3u8, 5, 7][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[2u8, 4, 6][..]).expect("valid polynomial");

        let (quot, rem) = a.div_rem(b).expect("division succeeds");

        // Quotient should be a constant (degree 0)
        assert_eq!(quot.degree(), 0);

        // Verify: quot * b + rem = a
        let product = quot.mul(&b).expect("multiplication succeeds");
        let reconstructed = (product + rem).expect("addition succeeds");

        assert_eq!(reconstructed, a);
    }

    #[test]
    fn div_rem_self_division() {
        // a / a = (1, 0)
        let a = Polynomial::try_from(&[3u8, 5, 7, 11][..]).expect("valid polynomial");

        let (quot, rem) = a.div_rem(a).expect("division succeeds");

        assert_eq!(quot.coefficients(), &[1]);
        assert_eq!(rem.coefficients(), &[0]);
    }

    #[test]
    fn div_rem_divisor_with_internal_zeros() {
        // Divisor like x^3 + 1 (coefficients [1, 0, 0, 1])
        let a = Polynomial::try_from(&[1u8, 2, 3, 4, 5][..]).expect("valid polynomial");
        let b: &[u8] = &[1, 0, 0, 1];

        let (quot, rem) = a.div_rem(b).expect("division succeeds");

        // Verify: quot * b + rem = a
        let product = quot.mul(b).expect("multiplication succeeds");
        let reconstructed = (product + rem).expect("addition succeeds");

        assert_eq!(reconstructed, a);
    }

    #[test]
    fn div_rem_high_degree() {
        // Test with higher degree polynomials
        let mut a_coeffs = [0u8; 100];
        let mut b_coeffs = [0u8; 30];

        for (i, c) in a_coeffs.iter_mut().enumerate() {
            *c = u8::try_from(i + 1).expect("index fits in u8");
        }
        for (i, c) in b_coeffs.iter_mut().enumerate() {
            *c = u8::try_from(i + 5).expect("index fits in u8");
        }

        let a = Polynomial::try_from(&a_coeffs[..]).expect("valid polynomial");
        let b = Polynomial::try_from(&b_coeffs[..]).expect("valid polynomial");

        let (quot, rem) = a.div_rem(b).expect("division succeeds");

        // Verify: quot * b + rem = a
        let product = quot.mul(&b).expect("multiplication succeeds");
        let reconstructed = (product + rem).expect("addition succeeds");

        assert_eq!(reconstructed, a);
    }
}
