use crate::Polynomial;

impl Polynomial {
    /// Returns `true` if this is the zero polynomial.
    ///
    /// A polynomial is zero if all its coefficients are zero. Due to internal
    /// normalization, this is equivalent to having degree 0 and a zero constant term.
    #[must_use]
    #[inline]
    pub const fn is_zero(&self) -> bool {
        self.degree == 0 && self.coefficients[0] == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_polynomial_is_zero() {
        let zero = Polynomial::default();

        assert!(zero.is_zero());
    }

    #[test]
    fn constant_nonzero_is_not_zero() {
        let mut p = Polynomial::default();
        p.set(0, 1);

        assert!(!p.is_zero());
    }

    #[test]
    fn polynomial_with_higher_degree_is_not_zero() {
        let mut p = Polynomial::default();
        p.set(2, 1);

        assert!(!p.is_zero());
    }

    #[test]
    fn collected_zeros_is_zero() {
        let p: Polynomial = [0u8; 4].into_iter().collect();

        assert!(p.is_zero());
    }

    #[test]
    fn collected_nonzero_is_not_zero() {
        let p: Polynomial = [1u8, 0, 0, 0].into_iter().collect();

        assert!(!p.is_zero());
    }

    const _: () = {
        // Compile-time verification that is_zero is const
        const ZERO: Polynomial = Polynomial {
            coefficients: [0; 255],
            degree: 0,
        };

        assert!(ZERO.is_zero());
    };
}
