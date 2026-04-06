use std::hash::{Hash, Hasher};

use crate::Polynomial;

impl Hash for Polynomial {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.degree.hash(state);
        self.coefficients[..=self.degree as usize].hash(state);
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    use crate::Polynomial;

    fn hash_of<T: Hash>(value: &T) -> u64 {
        let mut hasher = DefaultHasher::new();
        value.hash(&mut hasher);
        hasher.finish()
    }

    #[test]
    fn equal_polynomials_same_hash() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");

        assert_eq!(hash_of(&a), hash_of(&b));
    }

    #[test]
    fn trailing_zeros_same_hash() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2, 3, 0, 0][..]).expect("valid");

        assert_eq!(a, b);
        assert_eq!(hash_of(&a), hash_of(&b));
    }

    #[test]
    fn different_polynomials_likely_different_hash() {
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 2, 4][..]).expect("valid");

        // Not guaranteed, but extremely likely with a good hasher
        assert_ne!(hash_of(&a), hash_of(&b));
    }

    #[test]
    fn zero_polynomial_consistent_hash() {
        let a = Polynomial::default();
        let b = Polynomial::default();

        assert_eq!(hash_of(&a), hash_of(&b));
    }

    #[test]
    fn different_degrees_different_hash() {
        let a = Polynomial::try_from(&[1u8][..]).expect("valid");
        let b = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid");

        assert_ne!(hash_of(&a), hash_of(&b));
    }
}
