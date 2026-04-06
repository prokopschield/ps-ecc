use crate::{finite_field::mul, Polynomial};

impl Polynomial {
    /// Fused multiply-XOR: `self ^= a * b`
    ///
    /// Computes the product of polynomials `a` and `b`, `XORing` the result into `self`.
    /// This avoids allocating an intermediate polynomial for the product.
    ///
    /// In GF(2^8), XOR is equivalent to both addition and subtraction, so this
    /// computes `self = self + a * b` or equivalently `self = self - a * b`.
    #[inline]
    pub fn mul_xor_assign(&mut self, a: &Self, b: &Self) {
        let a_deg = a.degree as usize;
        let b_deg = b.degree as usize;

        // Early exit if either operand is zero polynomial
        if a.coefficients[a_deg] == 0 || b.coefficients[b_deg] == 0 {
            return;
        }

        let product_deg = a_deg + b_deg;

        // Accumulate a * b into self via XOR
        for i in 0..=a_deg {
            let a_coef = a.coefficients[i];

            if a_coef == 0 {
                continue;
            }

            for j in 0..=b_deg {
                let b_coef = b.coefficients[j];

                if b_coef == 0 {
                    continue;
                }

                self.coefficients[i + j] ^= mul(a_coef, b_coef);
            }
        }

        // Update degree: max of original self degree and product degree
        #[allow(clippy::cast_possible_truncation)]
        if product_deg > self.degree as usize {
            self.degree = product_deg as u8;
        }

        self.trim_degree();
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn mul_xor_assign_basic() {
        // dest = [1], a = [1, 1], b = [1, 1]
        // a * b = [1, 0, 1] (since 1*1 ^ 1*1 = 0 for middle term in GF(2))
        // dest ^= a * b => [1] ^ [1, 0, 1] = [0, 0, 1]
        let mut dest = Polynomial::try_from(&[1u8][..]).expect("valid polynomial");
        let a = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        dest.mul_xor_assign(&a, &b);

        assert_eq!(dest.coefficients(), &[0, 0, 1]);
    }

    #[test]
    fn mul_xor_assign_with_zero_a() {
        let mut dest = Polynomial::try_from(&[5u8, 3][..]).expect("valid polynomial");
        let a = Polynomial::default(); // zero polynomial
        let b = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");

        dest.mul_xor_assign(&a, &b);

        // dest unchanged since a is zero
        assert_eq!(dest.coefficients(), &[5, 3]);
    }

    #[test]
    fn mul_xor_assign_with_zero_b() {
        let mut dest = Polynomial::try_from(&[5u8, 3][..]).expect("valid polynomial");
        let a = Polynomial::try_from(&[1u8, 2, 3][..]).expect("valid polynomial");
        let b = Polynomial::default(); // zero polynomial

        dest.mul_xor_assign(&a, &b);

        // dest unchanged since b is zero
        assert_eq!(dest.coefficients(), &[5, 3]);
    }

    #[test]
    fn mul_xor_assign_accumulates() {
        // Start with non-zero dest and verify XOR accumulation
        let mut dest = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial");
        let a = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial");

        // a * b = [1, 0, 1]
        // dest ^= [1, 0, 1] => [0, 0, 0] = 0
        dest.mul_xor_assign(&a, &b);

        assert_eq!(dest.coefficients(), &[0]);
    }

    #[test]
    fn mul_xor_assign_extends_degree() {
        // dest has degree 0, product has higher degree
        let mut dest = Polynomial::try_from(&[5u8][..]).expect("valid polynomial");
        let a = Polynomial::try_from(&[1u8, 0, 1][..]).expect("valid polynomial"); // x^2 + 1
        let b = Polynomial::try_from(&[1u8, 1][..]).expect("valid polynomial"); // x + 1

        // a * b = x^3 + x^2 + x + 1 = [1, 1, 1, 1]
        // dest ^= [1, 1, 1, 1] => [5^1, 1, 1, 1] = [4, 1, 1, 1]
        dest.mul_xor_assign(&a, &b);

        assert_eq!(dest.coefficients(), &[4, 1, 1, 1]);
    }
}
