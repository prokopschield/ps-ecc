use std::ops::MulAssign;

use crate::{finite_field::mul, Polynomial};

impl MulAssign<u8> for Polynomial {
    /// Multiplies all coefficients by a scalar in GF(2^8).
    ///
    /// If `scalar` is zero, the result is the zero polynomial.
    /// Otherwise, the degree is preserved since GF(2^8) has no zero divisors.
    fn mul_assign(&mut self, scalar: u8) {
        if scalar == 0 {
            self.coefficients.fill(0);
            self.degree = 0;

            return;
        }

        for coef in &mut self.coefficients[..=self.degree as usize] {
            *coef = mul(*coef, scalar);
        }
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn mul_assign_by_one() {
        let mut p = Polynomial::try_from(&[5u8, 10, 15][..]).expect("valid polynomial");
        let original = p;

        p *= 1;

        assert_eq!(p, original);
    }

    #[test]
    fn mul_assign_by_zero() {
        let mut p = Polynomial::try_from(&[5u8, 10, 15][..]).expect("valid polynomial");

        p *= 0;

        assert_eq!(p.degree(), 0);
        assert_eq!(p.coefficients(), &[0]);
    }

    #[test]
    fn mul_assign_by_two() {
        // In GF(2^8), 2 * x is just a field multiplication
        let mut p = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        p *= 2;

        // GF(2^8) multiplication: 1*2=2, 2*2=4, 3*2=6
        assert_eq!(p.coefficients(), &[2, 4, 6]);
    }

    #[test]
    fn mul_assign_preserves_degree() {
        let mut p = Polynomial::try_from(&[1u8, 0, 0, 5][..]).expect("valid polynomial");

        assert_eq!(p.degree(), 3);

        p *= 7;

        assert_eq!(p.degree(), 3);
    }

    #[test]
    fn mul_assign_zero_polynomial() {
        let mut p = Polynomial::default();

        p *= 42;

        assert_eq!(p.degree(), 0);
        assert_eq!(p.coefficients(), &[0]);
    }

    #[test]
    fn mul_assign_by_inverse_roundtrip() {
        use crate::finite_field::inv;

        let original = Polynomial::try_from(&[3u8, 7, 11][..]).expect("valid polynomial");
        let mut p = original;

        let scale = 5u8;
        let scale_inv = inv(scale).expect("5 is invertible").get();

        p *= scale;
        p *= scale_inv;

        assert_eq!(p, original);
    }

    #[test]
    fn mul_assign_distributive() {
        // (a + b) * c = a*c + b*c
        let a = Polynomial::try_from(&[1u8, 2][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[3u8, 4][..]).expect("valid polynomial");
        let c = 5u8;

        let mut sum = (a + b).expect("addition succeeds");
        sum *= c;

        let mut a_scaled = a;
        a_scaled *= c;

        let mut b_scaled = b;
        b_scaled *= c;

        let sum_of_scaled = (a_scaled + b_scaled).expect("addition succeeds");

        assert_eq!(sum, sum_of_scaled);
    }

    #[test]
    fn mul_assign_associative_with_field() {
        // (p * a) * b = p * (a * b) where a, b are field elements
        use crate::finite_field::mul;

        let original = Polynomial::try_from(&[2u8, 4, 6][..]).expect("valid polynomial");
        let a = 3u8;
        let b = 7u8;

        let mut left = original;
        left *= a;
        left *= b;

        let mut right = original;
        right *= mul(a, b);

        assert_eq!(left, right);
    }

    #[test]
    fn mul_assign_high_degree() {
        let mut coeffs = [0u8; 255];
        coeffs[0] = 1;
        coeffs[254] = 1;

        let mut p = Polynomial::try_from(&coeffs[..]).expect("valid polynomial");

        assert_eq!(p.degree(), 254);

        p *= 3;

        assert_eq!(p.degree(), 254);
        assert_eq!(p[0u8], 3);
        assert_eq!(p[254u8], 3);
    }

    #[test]
    fn mul_assign_idempotent_with_one() {
        let mut p = Polynomial::try_from(&[7u8, 13, 19][..]).expect("valid polynomial");

        p *= 1;
        p *= 1;
        p *= 1;

        assert_eq!(p.coefficients(), &[7, 13, 19]);
    }
}
