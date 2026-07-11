use std::ops::Mul;

use crate::{error::PolynomialMulError, Polynomial};

use super::mul_coefficient_slices;

impl Mul<&Polynomial> for &Polynomial {
    type Output = Result<Polynomial, PolynomialMulError>;

    /// Multiplies two polynomials in GF(256).
    ///
    /// # Errors
    ///
    /// Returns `DegreeOverflow` if the result degree exceeds 254.
    fn mul(self, rhs: &Polynomial) -> Self::Output {
        mul_coefficient_slices(self.coefficients(), rhs.coefficients())
    }
}

impl Mul<&Self> for Polynomial {
    type Output = Result<Self, PolynomialMulError>;

    fn mul(self, rhs: &Self) -> Self::Output {
        (&self).mul(rhs)
    }
}

impl Mul<Self> for Polynomial {
    type Output = Result<Self, PolynomialMulError>;

    fn mul(self, rhs: Self) -> Self::Output {
        (&self).mul(&rhs)
    }
}

impl Mul<Polynomial> for &Polynomial {
    type Output = Result<Polynomial, PolynomialMulError>;

    fn mul(self, rhs: Polynomial) -> Self::Output {
        self.mul(&rhs)
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::ops::Mul;

    use crate::{error::PolynomialMulError, Polynomial};

    #[test]
    fn mul_two_polynomials() {
        // (1 + 2x) * (3 + 4x) = 3 + 4x + 6x + 8x^2 = 3 + (4^6)x + 8x^2
        // In GF(256): 4 ^ 6 = 2
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[3u8, 4][..]).expect("valid polynomial");

        let result = a.mul(&b).expect("degree fits");

        // a[0]*b[0] = 1*3 = 3
        // a[0]*b[1] + a[1]*b[0] = 1*4 + 2*3 = 4 ^ 6 = 2
        // a[1]*b[1] = 2*4 = 8
        assert_eq!(result.coefficients(), &[3, 2, 8]);
    }

    #[test]
    fn mul_by_zero() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let zero = Polynomial::default();

        let result = a.mul(&zero).expect("degree fits");

        assert_eq!(result.degree(), 0);
        assert_eq!(result.coefficients(), &[0]);
    }

    #[test]
    fn mul_by_one() {
        let a = Polynomial::try_from(&[5u8, 10, 15][..]).expect("valid polynomial");
        let one = Polynomial::try_from(&[1u8][..]).expect("valid polynomial");

        let result = a.mul(&one).expect("degree fits");

        assert_eq!(result.coefficients(), &[5, 10, 15]);
    }

    #[test]
    fn mul_by_x() {
        // Multiplying by x shifts coefficients
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let x = Polynomial::try_from(&[0u8, 1][..]).expect("valid polynomial");

        let result = a.mul(&x).expect("degree fits");

        assert_eq!(result.coefficients(), &[0, 1, 2, 3]);
    }

    #[test]
    fn mul_degree_overflow() {
        // Create two polynomials whose product exceeds max degree
        let mut coeffs = [0u8; 200];

        coeffs[199] = 1; // degree 199

        let a = Polynomial::try_from(&coeffs[..]).expect("valid polynomial");

        // 199 + 199 = 398 > 254
        let result = a.mul(&a);

        assert_eq!(result, Err(PolynomialMulError::DegreeOverflow));
    }

    #[test]
    fn mul_with_references() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[3u8, 4][..]).expect("valid polynomial");

        // All four combinations
        let r1 = (&a).mul(&b).expect("degree fits");
        let r2 = (&a).mul(b).expect("degree fits");
        let r3 = a.mul(&b).expect("degree fits");
        let r4 = a.mul(b).expect("degree fits");

        assert_eq!(r1, r2);
        assert_eq!(r2, r3);
        assert_eq!(r3, r4);
    }

    #[test]
    fn mul_max_degree_exactly() {
        // 127 + 127 = 254, which is exactly MAX_DEGREE
        let mut coeffs = [0u8; 128];

        coeffs[127] = 1;

        let a = Polynomial::try_from(&coeffs[..]).expect("valid polynomial");

        let result = a.mul(&a).expect("should succeed at boundary");

        assert_eq!(result.degree(), 254);
    }

    #[test]
    fn mul_is_commutative() {
        let a = Polynomial::try_from(&[1u8, 2, 3, 4][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[5u8, 6, 7][..]).expect("valid polynomial");

        let ab = (&a).mul(&b).expect("degree fits");
        let ba = (&b).mul(&a).expect("degree fits");

        assert_eq!(ab, ba);
    }
}
