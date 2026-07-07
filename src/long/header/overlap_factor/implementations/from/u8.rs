use crate::long::header::OverlapFactor;

impl From<u8> for OverlapFactor {
    #[inline]
    fn from(byte: u8) -> Self {
        Self::from_u8(byte)
    }
}

#[cfg(test)]
mod tests {
    use crate::long::OverlapFactor;

    #[test]
    fn test_from_u8_impl_matches_from_u8() {
        for byte in [0x00, 0x40, 0x80, 0xC0, 0x3F, 0x7F, 0xBF, 0xFF] {
            assert_eq!(OverlapFactor::from(byte), OverlapFactor::from_u8(byte));
        }
    }

    #[test]
    fn test_into_overlap_factor() {
        let factor: OverlapFactor = 0x80u8.into();

        assert_eq!(factor, OverlapFactor::Triple);
    }
}
