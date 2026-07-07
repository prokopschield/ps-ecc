use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Number of parity bytes per segment: `2 * parity`, since each parity symbol spans two bytes.
    #[inline]
    #[must_use]
    pub const fn parity_bytes(&self) -> u8 {
        self.parity() << 1
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};
    use crate::MAX_PARITY;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_parity_bytes_doubles_parity() -> Result<(), TestError> {
        let header = LongEccHeader::new(3, OverlapFactor::Simple, 88, 50, 0)?;

        assert_eq!(header.parity_bytes(), 6);

        Ok(())
    }

    #[test]
    fn test_parity_bytes_zero() -> Result<(), TestError> {
        let header = LongEccHeader::new(0, OverlapFactor::Simple, 82, 50, 0)?;

        assert_eq!(header.parity_bytes(), 0);

        Ok(())
    }

    #[test]
    fn test_parity_bytes_ignores_overlap_bits() -> Result<(), TestError> {
        let header = LongEccHeader::new(5, OverlapFactor::Quadruple, 92, 50, 0)?;

        assert_eq!(header.parity_bytes(), 10);

        Ok(())
    }

    #[test]
    fn test_parity_bytes_maximum() -> Result<(), TestError> {
        let header = LongEccHeader::new(MAX_PARITY, OverlapFactor::Simple, 208, 50, 0)?;

        assert_eq!(header.parity_bytes(), 126);

        Ok(())
    }
}
