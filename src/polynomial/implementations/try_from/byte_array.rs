use crate::{Polynomial, PolynomialFromSliceError};

impl<const N: usize> TryFrom<[u8; N]> for Polynomial {
    type Error = PolynomialFromSliceError;

    fn try_from(value: [u8; N]) -> Result<Self, Self::Error> {
        Self::from_slice(&value)
    }
}

impl<const N: usize> TryFrom<&[u8; N]> for Polynomial {
    type Error = PolynomialFromSliceError;

    fn try_from(value: &[u8; N]) -> Result<Self, Self::Error> {
        Self::from_slice(value)
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn from_array() {
        let poly = Polynomial::try_from([1, 2, 3]).expect("array should succeed");

        assert_eq!(poly.coefficients(), &[1, 2, 3]);
    }

    #[test]
    fn from_array_ref() {
        let poly = Polynomial::try_from(&[1, 2, 3]).expect("array ref should succeed");

        assert_eq!(poly.coefficients(), &[1, 2, 3]);
    }
}
