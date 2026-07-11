use crate::Polynomial;

pub const POLYNOMIAL_MAX_COEFFICIENTS_U8: u8 = u8::MAX;
pub const POLYNOMIAL_MAX_DEGREE_U8: u8 = POLYNOMIAL_MAX_COEFFICIENTS_U8 - 1;

pub const POLYNOMIAL_MAX_COEFFICIENTS: usize = POLYNOMIAL_MAX_COEFFICIENTS_U8 as usize;

impl Polynomial {
    /// Maximum degree of a polynomial over GF(256).
    pub const MAX_DEGREE: u8 = POLYNOMIAL_MAX_DEGREE_U8;

    /// Maximum number of coefficients in a polynomial.
    pub const MAX_COEFFICIENTS: u8 = POLYNOMIAL_MAX_COEFFICIENTS_U8;
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn max_degree_is_254() {
        assert_eq!(Polynomial::MAX_DEGREE, 254);
    }

    #[test]
    fn max_coefficients_is_255() {
        assert_eq!(Polynomial::MAX_COEFFICIENTS, 255);
    }

    #[test]
    fn max_coefficients_is_max_degree_plus_one() {
        assert_eq!(Polynomial::MAX_COEFFICIENTS, Polynomial::MAX_DEGREE + 1);
    }
}
