use crate::{error::PolynomialSetCoefficientsError, Polynomial};

impl Polynomial {
    /// Sets coefficients starting at degree `offset` from a slice.
    ///
    /// Copies `values` into the coefficient array at `offset..offset + values.len()`,
    /// then adjusts the degree accordingly.
    ///
    /// Empty slices are a no-op.
    ///
    /// # Errors
    ///
    /// Returns `OutOfBounds` if `offset + values.len()` exceeds 255.
    pub fn set_coefficients(
        &mut self,
        offset: u8,
        values: &[u8],
    ) -> Result<(), PolynomialSetCoefficientsError> {
        if values.is_empty() {
            return Ok(());
        }

        let start = offset as usize;
        let end = start + values.len();

        if end > Self::MAX_COEFFICIENTS as usize {
            return Err(PolynomialSetCoefficientsError::OutOfBounds { offset, end });
        }

        self.coefficients[start..end].copy_from_slice(values);

        #[allow(clippy::cast_possible_truncation)]
        let last = (end - 1) as u8;

        if last > self.degree {
            self.degree = last;
        }

        self.trim_degree();

        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::{error::PolynomialSetCoefficientsError, Polynomial};

    #[test]
    fn below_degree() {
        let mut p = Polynomial::try_from(&[1u8, 2, 3, 4, 5][..]).expect("valid polynomial");

        p.set_coefficients(1, &[10, 20]).expect("in bounds");

        assert_eq!(p.coefficients(), &[1, 10, 20, 4, 5]);
    }

    #[test]
    fn extends_degree() {
        let mut p = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");

        p.set_coefficients(2, &[3, 4, 5]).expect("in bounds");

        assert_eq!(p.coefficients(), &[1, 2, 3, 4, 5]);
    }

    #[test]
    fn trims_degree() {
        let mut p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        p.set_coefficients(1, &[0, 0]).expect("in bounds");

        assert_eq!(p.coefficients(), &[1]);
        assert_eq!(p.degree(), 0);
    }

    #[test]
    fn trims_to_zero_polynomial() {
        let mut p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        p.set_coefficients(0, &[0, 0, 0]).expect("in bounds");

        assert_eq!(p.coefficients(), &[0]);
        assert_eq!(p.degree(), 0);
    }

    #[test]
    fn empty_is_noop() {
        let mut p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let before = p;

        p.set_coefficients(0, &[]).expect("in bounds");

        assert_eq!(p, before);
    }

    #[test]
    fn at_offset_zero() {
        let mut p = Polynomial::try_from(&[0u8, 0, 5][..]).expect("valid polynomial");

        p.set_coefficients(0, &[1, 2]).expect("in bounds");

        assert_eq!(p.coefficients(), &[1, 2, 5]);
    }

    #[test]
    fn single_element() {
        let mut p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        p.set_coefficients(1, &[99]).expect("in bounds");

        assert_eq!(p.coefficients(), &[1, 99, 3]);
    }

    #[test]
    fn overwrites_entire_polynomial() {
        let mut p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        p.set_coefficients(0, &[10, 20, 30]).expect("in bounds");

        assert_eq!(p.coefficients(), &[10, 20, 30]);
    }

    #[test]
    fn gap_between_old_and_new() {
        let mut p = Polynomial::try_from(&[1u8][..]).expect("valid polynomial");

        p.set_coefficients(3, &[5, 6]).expect("in bounds");

        assert_eq!(p.coefficients(), &[1, 0, 0, 5, 6]);
    }

    #[test]
    fn extends_with_trailing_zeros() {
        let mut p = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");

        p.set_coefficients(3, &[5, 0]).expect("in bounds");

        assert_eq!(p.coefficients(), &[1, 2, 0, 5]);
        assert_eq!(p.degree(), 3);
    }

    #[test]
    fn overwrites_leading_coefficient() {
        let mut p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        p.set_coefficients(2, &[99]).expect("in bounds");

        assert_eq!(p.coefficients(), &[1, 2, 99]);
        assert_eq!(p.degree(), 2);
    }

    #[test]
    fn out_of_bounds_returns_error() {
        let mut p = Polynomial::default();

        let result = p.set_coefficients(250, &[1, 2, 3, 4, 5, 6]);

        assert_eq!(
            result,
            Err(PolynomialSetCoefficientsError::OutOfBounds {
                offset: 250,
                end: 256
            })
        );
    }
}
