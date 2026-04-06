use crate::Polynomial;

impl AsRef<[u8]> for Polynomial {
    fn as_ref(&self) -> &[u8] {
        self.coefficients()
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn as_ref_returns_coefficients() {
        let p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        let slice: &[u8] = p.as_ref();

        assert_eq!(slice, &[1, 2, 3]);
    }

    #[test]
    fn as_ref_zero_polynomial() {
        let p = Polynomial::default();

        let slice: &[u8] = p.as_ref();

        assert_eq!(slice, &[0]);
    }

    #[test]
    fn as_ref_returns_trimmed_slice() {
        // Input has trailing zeros, but coefficients() returns trimmed slice
        let p = Polynomial::try_from(&[1u8, 2, 3, 0, 0][..]).expect("valid polynomial");

        let slice: &[u8] = p.as_ref();

        assert_eq!(slice, &[1, 2, 3]);
    }

    #[test]
    fn as_ref_works_with_generic_fn() {
        fn takes_as_ref(x: impl AsRef<[u8]>) -> usize {
            x.as_ref().len()
        }

        let p = Polynomial::try_from(&[1u8, 2, 3, 4, 5][..]).expect("valid polynomial");

        assert_eq!(takes_as_ref(p), 5);
    }
}
