use std::fmt::{self, Debug};

use crate::Polynomial;

impl Debug for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Polynomial")
            .field("coefficients", &self.coefficients())
            .field("degree", &self.degree)
            .finish()
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn debug_shows_trimmed_coefficients() {
        let p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");

        let debug = format!("{p:?}");

        assert!(debug.contains("[1, 2, 3]"));
        assert!(!debug.contains("[1, 2, 3, 0"));
    }

    #[test]
    fn debug_shows_degree() {
        let p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");

        let debug = format!("{p:?}");

        assert!(debug.contains("degree: 2"));
    }

    #[test]
    fn debug_zero_polynomial() {
        let p = Polynomial::default();

        let debug = format!("{p:?}");

        assert!(debug.contains("[0]"));
        assert!(debug.contains("degree: 0"));
    }

    #[test]
    fn debug_alternate_format() {
        let p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");

        let debug = format!("{p:#?}");

        assert!(debug.contains("Polynomial"));
        assert!(debug.contains("coefficients"));
        assert!(debug.contains("degree"));
    }
}
