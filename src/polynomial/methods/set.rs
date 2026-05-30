use std::cmp::Ordering::{Equal, Greater, Less};

use crate::Polynomial;

impl Polynomial {
    /// Sets the coefficient of degree `idx` to `value`.
    ///
    /// If `value` is non-zero and `idx` exceeds the current degree, the degree
    /// increases to `idx`. If `value` is zero and `idx` equals the current degree,
    /// the degree decreases to the next non-zero coefficient.
    ///
    /// # Panics
    ///
    /// Panics if `idx` exceeds 254, the maximum degree of a polynomial in GF(256).
    pub fn set(&mut self, idx: u8, value: u8) {
        self.coefficients[idx as usize] = value;

        match (idx.cmp(&self.degree), value) {
            (Equal, 0) => self.trim_degree(),
            (Greater, 1..) => self.degree = idx,
            (Less, ..) | (Greater, 0) | (Equal, 1..) => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn set_below_degree() {
        let mut p = Polynomial {
            degree: 3,
            ..Default::default()
        };

        p.coefficients[3] = 5;

        p.set(1, 7);

        assert_eq!(p.coefficients[1], 7);
        assert_eq!(p.degree, 3);
    }

    #[test]
    fn set_at_degree_nonzero() {
        let mut p = Polynomial {
            degree: 2,
            ..Default::default()
        };

        p.coefficients[2] = 3;

        p.set(2, 9);

        assert_eq!(p.coefficients[2], 9);
        assert_eq!(p.degree, 2);
    }

    #[test]
    fn set_at_degree_zero_trims() {
        let mut p = Polynomial {
            degree: 3,
            ..Default::default()
        };

        p.coefficients[0] = 1;
        p.coefficients[1] = 2;
        p.coefficients[3] = 4;

        p.set(3, 0);

        assert_eq!(p.coefficients[3], 0);
        assert_eq!(p.degree, 1);
    }

    #[test]
    fn set_above_degree_nonzero() {
        let mut p = Polynomial {
            degree: 1,
            ..Default::default()
        };

        p.coefficients[1] = 1;

        p.set(5, 3);

        assert_eq!(p.coefficients[5], 3);
        assert_eq!(p.degree, 5);
    }

    #[test]
    fn set_above_degree_zero() {
        let mut p = Polynomial {
            degree: 2,
            ..Default::default()
        };

        p.coefficients[2] = 1;

        p.set(5, 0);

        assert_eq!(p.coefficients[5], 0);
        assert_eq!(p.degree, 2);
    }

    #[test]
    fn set_trims_to_zero_polynomial() {
        let mut p = Polynomial {
            degree: 0,
            ..Default::default()
        };

        p.coefficients[0] = 5;

        p.set(0, 0);

        assert_eq!(p.coefficients[0], 0);
        assert_eq!(p.degree, 0);
    }
}
