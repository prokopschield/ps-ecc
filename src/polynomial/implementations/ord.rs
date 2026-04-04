use std::cmp::Ordering;

use crate::Polynomial;

impl Ord for Polynomial {
    fn cmp(&self, other: &Self) -> Ordering {
        self.degree().cmp(&other.degree()).then_with(|| {
            let d = self.degree() as usize;
            self.coefficients[..=d]
                .iter()
                .rev()
                .cmp(other.coefficients[..=d].iter().rev())
        })
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering::{Equal, Greater, Less};

    use crate::Polynomial;

    #[test]
    fn equal_zero_polynomials() {
        let a = Polynomial::default();
        let b = Polynomial::default();
        assert_eq!(a.cmp(&b), Equal);
    }

    #[test]
    fn equal_nonzero_polynomials() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();
        a.set(2, 5);
        a.set(0, 3);
        b.set(2, 5);
        b.set(0, 3);
        assert_eq!(a.cmp(&b), Equal);
    }

    #[test]
    fn higher_degree_is_greater() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();
        a.set(3, 1);
        b.set(2, 255);
        assert_eq!(a.cmp(&b), Greater);
        assert_eq!(b.cmp(&a), Less);
    }

    #[test]
    fn same_degree_compare_leading_coefficient() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();
        a.set(2, 10);
        b.set(2, 5);
        assert_eq!(a.cmp(&b), Greater);
        assert_eq!(b.cmp(&a), Less);
    }

    #[test]
    fn same_degree_equal_leading_compare_next() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();
        a.set(2, 5);
        a.set(1, 10);
        b.set(2, 5);
        b.set(1, 3);
        assert_eq!(a.cmp(&b), Greater);
        assert_eq!(b.cmp(&a), Less);
    }

    #[test]
    fn same_degree_differ_only_at_constant() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();
        a.set(2, 5);
        a.set(0, 2);
        b.set(2, 5);
        b.set(0, 1);
        assert_eq!(a.cmp(&b), Greater);
        assert_eq!(b.cmp(&a), Less);
    }

    #[test]
    fn reflexive() {
        let mut p = Polynomial::default();
        p.set(3, 7);
        p.set(1, 2);
        assert_eq!(p.cmp(&p), Equal);
    }

    #[test]
    fn transitive() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();
        let mut c = Polynomial::default();
        a.set(1, 1);
        b.set(1, 2);
        c.set(1, 3);
        assert_eq!(a.cmp(&b), Less);
        assert_eq!(b.cmp(&c), Less);
        assert_eq!(a.cmp(&c), Less);
    }
}
