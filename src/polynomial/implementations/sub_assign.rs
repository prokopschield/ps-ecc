use std::ops::{BitXorAssign, SubAssign};

use crate::Polynomial;

impl SubAssign for Polynomial {
    fn sub_assign(&mut self, rhs: Self) {
        BitXorAssign::bitxor_assign(self, rhs);
    }
}

impl SubAssign<&Self> for Polynomial {
    fn sub_assign(&mut self, rhs: &Self) {
        BitXorAssign::bitxor_assign(self, rhs);
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn sub_assign_basic() {
        let mut a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        a -= b;

        assert_eq!(a.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn sub_assign_different_degrees() {
        let mut a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[3u8, 4, 5, 6][..]).expect("valid polynomial");

        a -= b;

        assert_eq!(a.coefficients(), &[1 ^ 3, 2 ^ 4, 5, 6]);
    }

    #[test]
    fn sub_assign_cancels_to_zero() {
        let mut a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        a -= b;

        assert_eq!(a.degree(), 0);
        assert_eq!(a.coefficients(), &[0]);
    }

    #[test]
    fn sub_assign_with_reference() {
        let mut a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        a -= &b;

        assert_eq!(a.coefficients(), &[1 ^ 4, 2 ^ 5, 3 ^ 6]);
    }

    #[test]
    fn sub_assign_multiple() {
        let mut a = Polynomial::try_from(&[1u8][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[2u8][..]).expect("valid polynomial");
        let c = Polynomial::try_from(&[4u8][..]).expect("valid polynomial");

        a -= b;
        a -= c;

        assert_eq!(a.coefficients(), &[1 ^ 2 ^ 4]);
    }

    #[test]
    fn sub_assign_equals_add_assign() {
        let mut a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let mut b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let c = Polynomial::try_from(&[4u8, 5, 6][..]).expect("valid polynomial");

        a -= &c;
        b += &c;

        assert_eq!(a, b);
    }
}
