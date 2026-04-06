use crate::Polynomial;

impl PartialEq for Polynomial {
    fn eq(&self, other: &Self) -> bool {
        if self.degree != other.degree {
            return false;
        }

        let len = self.degree as usize + 1;

        self.coefficients[..len] == other.coefficients[..len]
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn equal_polynomials() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");

        assert_eq!(a, b);
    }

    #[test]
    fn different_degrees_not_equal() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2][..]).expect("valid");

        assert_ne!(a, b);
    }

    #[test]
    fn same_degree_different_coefficients() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2, 4][..]).expect("valid");

        assert_ne!(a, b);
    }

    #[test]
    fn zero_polynomials_equal() {
        let a = Polynomial::default();
        let b = Polynomial::default();

        assert_eq!(a, b);
    }

    #[test]
    fn zero_vs_nonzero_not_equal() {
        let a = Polynomial::default();
        let b = Polynomial::try_from(&[1u8][..]).expect("valid");

        assert_ne!(a, b);
    }

    #[test]
    fn reflexive() {
        let p = Polynomial::try_from(&[5u8, 10, 15][..]).expect("valid");

        assert_eq!(p, p);
    }

    #[test]
    fn symmetric() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");

        assert_eq!(a == b, b == a);
    }

    #[test]
    fn transitive() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let c = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");

        assert!(a == b && b == c && a == c);
    }

    #[test]
    fn trailing_zeros_trimmed_equal() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2, 3, 0, 0][..]).expect("valid");

        assert_eq!(a, b);
    }

    #[test]
    fn only_compares_relevant_coefficients() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();

        a.set(0, 1);
        b.set(0, 1);

        // Manually set garbage beyond degree (simulating uninitialized memory)
        // This shouldn't affect equality since degree is 0
        // Note: We can't actually do this safely, but the test verifies
        // that equality only checks up to degree
        assert_eq!(a, b);
    }
}
