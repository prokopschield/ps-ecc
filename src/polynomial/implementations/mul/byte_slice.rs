use std::ops::Mul;

use crate::{error::PolynomialMulError, Polynomial};

use super::mul_coefficient_slices;

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

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::ops::Mul;

    use crate::{error::PolynomialMulError, Polynomial};

    #[test]
    fn mul_by_slice() {
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b: &[u8] = &[3, 4];

        let result = a.mul(b).expect("degree fits");

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
