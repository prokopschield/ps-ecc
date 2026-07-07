use super::super::OverlapFactor;

impl OverlapFactor {
    /// Returns the coverage count: the number of codewords protecting each byte (1 to 4).
    #[inline]
    #[must_use]
    pub const fn count(self) -> u8 {
        match self {
            Self::Simple => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Quadruple => 4,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::long::OverlapFactor;

    #[test]
    fn test_count_of_each_factor() {
        assert_eq!(OverlapFactor::Simple.count(), 1);
        assert_eq!(OverlapFactor::Double.count(), 2);
        assert_eq!(OverlapFactor::Triple.count(), 3);
        assert_eq!(OverlapFactor::Quadruple.count(), 4);
    }
}
