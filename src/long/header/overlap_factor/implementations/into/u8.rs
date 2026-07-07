use crate::long::header::OverlapFactor;

impl From<OverlapFactor> for u8 {
    #[inline]
    fn from(factor: OverlapFactor) -> Self {
        factor.to_u8()
    }
}

#[cfg(test)]
mod tests {
    use crate::long::OverlapFactor;

    #[test]
    fn test_u8_from_overlap_factor() {
        assert_eq!(u8::from(OverlapFactor::Simple), 0);
        assert_eq!(u8::from(OverlapFactor::Double), 64);
        assert_eq!(u8::from(OverlapFactor::Triple), 128);
        assert_eq!(u8::from(OverlapFactor::Quadruple), 192);
    }

    #[test]
    fn test_conversion_roundtrip() {
        for factor in [
            OverlapFactor::Simple,
            OverlapFactor::Double,
            OverlapFactor::Triple,
            OverlapFactor::Quadruple,
        ] {
            assert_eq!(OverlapFactor::from(u8::from(factor)), factor);
        }
    }
}
