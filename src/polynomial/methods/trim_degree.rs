use crate::{polynomial::constants::POLYNOMIAL_MAX_DEGREE_U8, Polynomial};

impl Polynomial {
    pub(crate) const fn trim_degree(&mut self) {
        if self.degree > POLYNOMIAL_MAX_DEGREE_U8 {
            self.degree = POLYNOMIAL_MAX_DEGREE_U8;
        }

        while self.degree > 0 && self.coefficients[self.degree as usize] == 0 {
            self.degree -= 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn trim_single_trailing_zero() {
        let mut p = Polynomial {
            degree: 3,
            ..Default::default()
        };
        p.coefficients[0] = 1;
        p.coefficients[1] = 2;
        p.coefficients[2] = 3;
        // degree 3 coefficient is 0

        p.trim_degree();
        assert_eq!(p.degree, 2);
    }

    #[test]
    fn trim_multiple_trailing_zeros() {
        let mut p = Polynomial {
            degree: 5,
            ..Default::default()
        };
        p.coefficients[0] = 1;
        // coefficients 1-5 are all 0

        p.trim_degree();
        assert_eq!(p.degree, 0);
    }

    #[test]
    fn trim_no_trailing_zeros() {
        let mut p = Polynomial {
            degree: 2,
            ..Default::default()
        };
        p.coefficients[0] = 1;
        p.coefficients[1] = 2;
        p.coefficients[2] = 3;

        p.trim_degree();
        assert_eq!(p.degree, 2);
    }

    #[test]
    fn trim_zero_polynomial() {
        let mut p = Polynomial::default();

        p.trim_degree();
        assert_eq!(p.degree, 0);
    }

    #[test]
    fn trim_all_zeros_except_constant() {
        let mut p = Polynomial {
            degree: 4,
            ..Default::default()
        };
        // all coefficients are 0

        p.trim_degree();
        assert_eq!(p.degree, 0);
    }
}
