use std::ops::Mul;

use crate::{
    error::PolynomialMulError,
    finite_field::{add, mul},
    Polynomial,
};

/// Core multiplication on trimmed coefficient slices.
///
/// Both slices must represent trimmed polynomials (no trailing zeros beyond degree).
/// Empty slices are treated as the zero polynomial.
fn mul_coefficient_slices(lhs: &[u8], rhs: &[u8]) -> Result<Polynomial, PolynomialMulError> {
    if lhs.is_empty() || rhs.is_empty() {
        return Ok(Polynomial::default());
    }

    let lhs_degree = lhs.len() - 1;
    let rhs_degree = rhs.len() - 1;

    let result_degree = lhs_degree + rhs_degree;

    if result_degree > Polynomial::MAX_DEGREE as usize {
        return Err(PolynomialMulError::DegreeOverflow);
    }

    #[allow(clippy::cast_possible_truncation)]
    let mut result = Polynomial {
        degree: result_degree as u8,
        ..Polynomial::default()
    };

    for (i, &lhs_coef) in lhs.iter().enumerate() {
        for (j, &rhs_coef) in rhs.iter().enumerate() {
            let k = i + j;

            result.coefficients[k] = add(result.coefficients[k], mul(lhs_coef, rhs_coef));
        }
    }

    result.trim_degree();

    Ok(result)
}

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

impl Mul<&[u8]> for &Polynomial {
    type Output = Result<Polynomial, PolynomialMulError>;

    /// Multiplies a polynomial by a coefficient slice in GF(256).
    ///
    /// # Errors
    ///
    /// Returns `DegreeOverflow` if the result degree exceeds 254.
    fn mul(self, rhs: &[u8]) -> Self::Output {
        let rhs_trimmed = match rhs.iter().rposition(|&c| c != 0) {
            Some(degree) => &rhs[..=degree],
            None => &[],
        };

        mul_coefficient_slices(self.coefficients(), rhs_trimmed)
    }
}

impl Mul<&[u8]> for Polynomial {
    type Output = Result<Self, PolynomialMulError>;

    fn mul(self, rhs: &[u8]) -> Self::Output {
        (&self).mul(rhs)
    }
}

impl<const N: usize> Mul<[u8; N]> for &Polynomial {
    type Output = Result<Polynomial, PolynomialMulError>;

    fn mul(self, rhs: [u8; N]) -> Self::Output {
        self.mul(rhs.as_slice())
    }
}

impl<const N: usize> Mul<[u8; N]> for Polynomial {
    type Output = Result<Self, PolynomialMulError>;

    fn mul(self, rhs: [u8; N]) -> Self::Output {
        (&self).mul(rhs.as_slice())
    }
}

impl<const N: usize> Mul<&[u8; N]> for &Polynomial {
    type Output = Result<Polynomial, PolynomialMulError>;

    fn mul(self, rhs: &[u8; N]) -> Self::Output {
        self.mul(rhs.as_slice())
    }
}

impl<const N: usize> Mul<&[u8; N]> for Polynomial {
    type Output = Result<Self, PolynomialMulError>;

    fn mul(self, rhs: &[u8; N]) -> Self::Output {
        (&self).mul(rhs.as_slice())
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

    #[test]
    fn mul_by_slice() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b: &[u8] = &[3, 4];

        let result = a.mul(b).expect("degree fits");

        assert_eq!(result.coefficients(), &[3, 2, 8]);
    }

    #[test]
    fn mul_by_array() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");

        let result = a.mul([3u8, 4]).expect("degree fits");

        assert_eq!(result.coefficients(), &[3, 2, 8]);
    }

    #[test]
    fn mul_by_array_ref() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");

        let result = a.mul(&[3u8, 4]).expect("degree fits");

        assert_eq!(result.coefficients(), &[3, 2, 8]);
    }

    #[test]
    fn mul_slice_degree_overflow() {
        let a = Polynomial::try_from(&[1u8][..]).expect("valid polynomial");
        let mut long_slice = [0u8; 256];

        long_slice[255] = 1; // degree 255, which exceeds MAX_DEGREE

        let result = a.mul(long_slice.as_slice());

        assert_eq!(result, Err(PolynomialMulError::DegreeOverflow));
    }

    #[test]
    fn mul_long_slice_with_trailing_zeros() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let mut long_slice = [0u8; 300];

        long_slice[0] = 3;
        long_slice[1] = 4; // effective degree is 1, despite slice length 300

        let result = a.mul(long_slice.as_slice()).expect("degree fits");

        assert_eq!(result.coefficients(), &[3, 2, 8]);
    }

    #[test]
    fn mul_slice_matches_polynomial() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        let via_poly = (&a).mul(&b).expect("degree fits");
        let via_slice = (&a).mul(b.coefficients()).expect("degree fits");

        assert_eq!(via_poly, via_slice);
    }
}
