use std::ops::Mul;

use crate::{error::PolynomialMulError, Polynomial};

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

    use crate::Polynomial;

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
}
