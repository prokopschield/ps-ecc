use super::super::OverlapFactor;

impl OverlapFactor {
    /// Extracts the overlap factor from the high two bits of a parity byte.
    ///
    /// The low six bits are ignored, so any byte maps to a valid factor.
    #[inline]
    #[must_use]
    pub const fn from_u8(byte: u8) -> Self {
        match byte & 0xC0 {
            0x00 => Self::Simple,
            0x40 => Self::Double,
            0x80 => Self::Triple,
            0xC0 => Self::Quadruple,
            // `& 0xC0` yields only the four values above
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::long::OverlapFactor;

    #[test]
    fn test_from_u8_high_bits() {
        assert_eq!(OverlapFactor::from_u8(0x00), OverlapFactor::Simple);
        assert_eq!(OverlapFactor::from_u8(0x40), OverlapFactor::Double);
        assert_eq!(OverlapFactor::from_u8(0x80), OverlapFactor::Triple);
        assert_eq!(OverlapFactor::from_u8(0xC0), OverlapFactor::Quadruple);
    }

    #[test]
    fn test_from_u8_ignores_low_bits() {
        assert_eq!(OverlapFactor::from_u8(0x3F), OverlapFactor::Simple);
        assert_eq!(OverlapFactor::from_u8(0x7F), OverlapFactor::Double);
        assert_eq!(OverlapFactor::from_u8(0xBF), OverlapFactor::Triple);
        assert_eq!(OverlapFactor::from_u8(0xFF), OverlapFactor::Quadruple);
    }
}
