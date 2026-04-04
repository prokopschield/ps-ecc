use crate::{Polynomial, PolynomialFromSliceError};

impl Polynomial {
    /// Creates a polynomial from a slice of coefficients.
    ///
    /// Coefficients are ordered from degree 0 (constant term) to the leading term.
    ///
    /// # Errors
    ///
    /// Returns [`PolynomialFromSliceError::TooLong`] if the slice exceeds 255 elements.
    pub fn from_slice(coeffs: &[u8]) -> Result<Self, PolynomialFromSliceError> {
        let size = coeffs.len();
        let size_u8 = u8::try_from(size).map_err(|_| PolynomialFromSliceError::TooLong { size })?;
        let degree = size_u8.saturating_sub(1);

        let mut poly = Self::default();

        poly.coefficients[..size].copy_from_slice(coeffs);
        poly.degree = degree;
        poly.trim_degree();

        Ok(poly)
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::{Polynomial, PolynomialFromSliceError};

    #[test]
    fn empty_slice() {
        let poly = Polynomial::from_slice(&[]).expect("empty slice should succeed");

        assert_eq!(poly.degree(), 0);
        assert_eq!(poly.coefficients(), &[0]);
    }

    #[test]
    fn single_coefficient() {
        let poly = Polynomial::from_slice(&[42]).expect("single coefficient should succeed");

        assert_eq!(poly.degree(), 0);
        assert_eq!(poly.coefficients(), &[42]);
    }

    #[test]
    fn multiple_coefficients() {
        let poly =
            Polynomial::from_slice(&[1, 2, 3]).expect("multiple coefficients should succeed");

        assert_eq!(poly.degree(), 2);
        assert_eq!(poly.coefficients(), &[1, 2, 3]);
    }

    #[test]
    fn trims_leading_zeros() {
        let poly =
            Polynomial::from_slice(&[1, 2, 0, 0]).expect("slice with leading zeros should succeed");

        assert_eq!(poly.degree(), 1);
        assert_eq!(poly.coefficients(), &[1, 2]);
    }

    #[test]
    fn all_zeros() {
        let poly = Polynomial::from_slice(&[0, 0, 0]).expect("all zeros should succeed");

        assert_eq!(poly.degree(), 0);
        assert_eq!(poly.coefficients(), &[0]);
    }

    #[test]
    fn sparse_polynomial() {
        let poly = Polynomial::from_slice(&[1, 0, 0, 5]).expect("sparse polynomial should succeed");

        assert_eq!(poly.degree(), 3);
        assert_eq!(poly.coefficients(), &[1, 0, 0, 5]);
    }

    #[test]
    fn max_size() {
        let coeffs = [1u8; 255];
        let poly = Polynomial::from_slice(&coeffs).expect("max size should succeed");

        assert_eq!(poly.degree(), 254);
        assert_eq!(poly.coefficients().len(), 255);
    }

    #[test]
    fn too_long() {
        let coeffs = [1u8; 256];
        let err = Polynomial::from_slice(&coeffs).expect_err("256 coefficients should fail");

        assert!(matches!(
            err,
            PolynomialFromSliceError::TooLong { size: 256 }
        ));
    }

    #[test]
    fn roundtrip_with_coefficients() {
        let original = &[5, 0, 3, 7];
        let poly = Polynomial::from_slice(original).expect("roundtrip should succeed");

        assert_eq!(poly.coefficients(), original);
    }
}
