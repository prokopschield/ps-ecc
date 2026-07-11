use crate::{Polynomial, PolynomialFromSliceError};

impl TryFrom<&[u8]> for Polynomial {
    type Error = PolynomialFromSliceError;

    fn try_from(value: &[u8]) -> Result<Self, Self::Error> {
        Self::from_slice(value)
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn from_slice() {
        let slice: &[u8] = &[1, 2, 3];
        let poly = Polynomial::try_from(slice).expect("slice should succeed");

        assert_eq!(poly.coefficients(), &[1, 2, 3]);
    }

    #[test]
    fn from_vec_ref() {
        let vec = vec![1, 2, 3];
        let poly = Polynomial::try_from(vec.as_slice()).expect("vec ref should succeed");

        assert_eq!(poly.coefficients(), &[1, 2, 3]);
    }
}
