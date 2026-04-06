use crate::{
    finite_field::{inv, mul, sub},
    Polynomial, PolynomialDivError,
};

impl Polynomial {
    /// In-place polynomial division: `self` becomes the remainder, quotient is output.
    ///
    /// After this method, `self` contains the remainder and `quotient` contains
    /// the quotient of dividing the original `self` by `divisor`.
    ///
    /// # Errors
    ///
    /// Returns `PolynomialDivError::ZeroDivisor` if the divisor is the zero polynomial.
    /// Returns `PolynomialDivError::GFError` if field inversion fails.
    #[inline]
    pub fn div_rem_inplace(
        &mut self,
        divisor: &Self,
        quotient: &mut Self,
    ) -> Result<(), PolynomialDivError> {
        let divisor_deg = divisor.degree as usize;

        // Check for zero divisor
        if divisor_deg == 0 && divisor.coefficients[0] == 0 {
            return Err(PolynomialDivError::ZeroDivisor);
        }

        // Precompute inverse of leading coefficient
        let divisor_lc_inv = inv(divisor.coefficients[divisor_deg])?.get();

        // Clear quotient
        quotient.coefficients.fill(0);
        quotient.degree = 0;

        let mut deg = self.degree as usize;

        // If dividend degree < divisor degree, remainder is dividend, quotient is 0
        if deg < divisor_deg {
            return Ok(());
        }

        // Set quotient degree
        #[allow(clippy::cast_possible_truncation)]
        {
            quotient.degree = (deg - divisor_deg) as u8;
        }

        while deg >= divisor_deg {
            let lead_coef = self.coefficients[deg];

            if lead_coef == 0 {
                if deg == 0 {
                    break;
                }

                deg -= 1;
                continue;
            }

            let ratio = mul(lead_coef, divisor_lc_inv);
            let shift = deg - divisor_deg;

            quotient.coefficients[shift] = ratio;
            self.coefficients[deg] = 0;

            // Subtract ratio * divisor from self (in GF(2^8), sub == xor)
            for i in 0..divisor_deg {
                self.coefficients[shift + i] = sub(
                    self.coefficients[shift + i],
                    mul(ratio, divisor.coefficients[i]),
                );
            }

            if deg == 0 {
                break;
            }

            deg -= 1;
        }

        // Update degrees
        self.trim_degree();
        quotient.trim_degree();

        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::{finite_field::mul, Polynomial};

    #[test]
    fn div_rem_inplace_basic() {
        // (x^2 + 1) / (x + 1) in GF(256)
        // = (x + 1) with remainder 0
        let mut dividend = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");

        let divisor = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");
        let mut quotient = Polynomial::default();

        dividend
            .div_rem_inplace(&divisor, &mut quotient)
            .expect("division succeeds");

        assert_eq!(quotient.coefficients(), &[1, 1]);
        assert_eq!(dividend.coefficients(), &[0]);
    }

    #[test]
    fn div_rem_inplace_with_remainder() {
        // (x^2 + x + 1) / (x + 1) should have a remainder
        let mut dividend = Polynomial::try_from(&[1u8, 1, 1][..]).expect("valid polynomial");

        let divisor = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        let original = Polynomial::try_from(&[1u8, 1, 1][..]).expect("valid polynomial");

        let mut quotient = Polynomial::default();

        dividend
            .div_rem_inplace(&divisor, &mut quotient)
            .expect("division succeeds");

        // Verify: quotient * divisor + remainder = original
        let mut reconstructed = Polynomial::default();

        // Multiply quotient * divisor
        for i in 0..=quotient.degree as usize {
            for j in 0..=divisor.degree as usize {
                reconstructed.coefficients[i + j] ^=
                    mul(quotient.coefficients[i], divisor.coefficients[j]);
            }
        }

        reconstructed.degree = quotient.degree + divisor.degree;
        reconstructed.trim_degree();

        // Add remainder
        for i in 0..=dividend.degree as usize {
            reconstructed.coefficients[i] ^= dividend.coefficients[i];
        }

        reconstructed.trim_degree();

        assert_eq!(reconstructed, original);
    }

    #[test]
    fn div_rem_inplace_smaller_dividend() {
        // dividend degree < divisor degree => remainder = dividend, quotient = 0
        let mut dividend = Polynomial::try_from(&[5u8, 3][..]).expect("valid polynomial");

        let divisor = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let original = Polynomial::try_from(&[5u8, 3][..]).expect("valid polynomial");
        let mut quotient = Polynomial::default();

        dividend
            .div_rem_inplace(&divisor, &mut quotient)
            .expect("division succeeds");

        assert_eq!(quotient.coefficients(), &[0]);
        assert_eq!(dividend, original);
    }

    #[test]
    fn div_rem_inplace_zero_divisor() {
        let mut dividend = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let divisor = Polynomial::default();
        let mut quotient = Polynomial::default();

        let result = dividend.div_rem_inplace(&divisor, &mut quotient);

        assert!(result.is_err());
    }
}
