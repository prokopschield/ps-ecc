use super::super::OverlapFactor;

impl OverlapFactor {
    /// Returns the bit-packed representation, occupying the high two bits of the parity byte.
    #[inline]
    #[must_use]
    pub const fn to_u8(self) -> u8 {
        self as u8
    }
}

#[cfg(test)]
mod tests {
    use crate::long::OverlapFactor;

    #[test]
    fn test_to_u8_packs_high_bits() {
        assert_eq!(OverlapFactor::Simple.to_u8(), 0);
        assert_eq!(OverlapFactor::Double.to_u8(), 64);
        assert_eq!(OverlapFactor::Triple.to_u8(), 128);
        assert_eq!(OverlapFactor::Quadruple.to_u8(), 192);
    }

    #[test]
    fn test_to_u8_from_u8_roundtrip() {
        for factor in [
            OverlapFactor::Simple,
            OverlapFactor::Double,
            OverlapFactor::Triple,
            OverlapFactor::Quadruple,
        ] {
            assert_eq!(OverlapFactor::from_u8(factor.to_u8()), factor);
        }
    }
}
