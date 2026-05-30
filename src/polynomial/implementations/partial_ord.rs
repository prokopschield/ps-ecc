use crate::Polynomial;

impl PartialOrd for Polynomial {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use crate::Polynomial;

    #[test]
    fn partial_cmp_returns_some() {
        let a = Polynomial::default();
        let b = Polynomial::default();

        assert!(a.partial_cmp(&b).is_some());
    }

    #[test]
    fn partial_cmp_consistent_with_cmp() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();

        a.set(2, 5);
        b.set(2, 3);

        assert_eq!(a.partial_cmp(&b), Some(a.cmp(&b)));
    }

    #[test]
    fn less_than_operator() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();

        a.set(1, 1);
        b.set(1, 2);

        assert!(a < b);
        assert!((b >= a));
    }

    #[test]
    fn greater_than_operator() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();

        a.set(2, 1);
        b.set(1, 255);

        assert!(a > b);
        assert!((b <= a));
    }

    #[test]
    fn less_or_equal_operator() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();

        a.set(1, 5);
        b.set(1, 5);

        assert!(a <= b);
        assert!(b <= a);

        b.set(1, 6);

        assert!(a <= b);
        assert!((b > a));
    }

    #[test]
    fn greater_or_equal_operator() {
        let mut a = Polynomial::default();
        let mut b = Polynomial::default();

        a.set(1, 5);
        b.set(1, 5);

        assert!(a >= b);
        assert!(b >= a);

        a.set(1, 6);

        assert!(a >= b);
        assert!((b < a));
    }
}
