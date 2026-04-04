use std::ops::Index;

use crate::Polynomial;

impl Index<usize> for Polynomial {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        self.coefficients.get(index).unwrap_or(&0)
    }
}

impl Index<u8> for Polynomial {
    type Output = u8;

    fn index(&self, index: u8) -> &Self::Output {
        &self[index as usize]
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn index_usize_within_bounds() {
        let p = Polynomial::from_slice(&[3, 7, 0, 5]).expect("valid polynomial");

        assert_eq!(p[0usize], 3);
        assert_eq!(p[1usize], 7);
        assert_eq!(p[2usize], 0);
        assert_eq!(p[3usize], 5);
    }

    #[test]
    fn index_usize_out_of_bounds_returns_zero() {
        let p = Polynomial::from_slice(&[1, 2]).expect("valid polynomial");

        assert_eq!(p[2usize], 0);
        assert_eq!(p[100usize], 0);
    }

    #[test]
    fn index_usize_empty_polynomial() {
        let p = Polynomial::from_slice(&[]).expect("valid polynomial");

        assert_eq!(p[0usize], 0);
    }

    #[test]
    fn index_u8_within_bounds() {
        let p = Polynomial::from_slice(&[10, 20, 30]).expect("valid polynomial");

        assert_eq!(p[0u8], 10);
        assert_eq!(p[1u8], 20);
        assert_eq!(p[2u8], 30);
    }

    #[test]
    fn index_u8_out_of_bounds_returns_zero() {
        let p = Polynomial::from_slice(&[1]).expect("valid polynomial");

        assert_eq!(p[1u8], 0);
        assert_eq!(p[u8::MAX], 0);
    }
}
