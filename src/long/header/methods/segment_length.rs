use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Usable data bytes per 255-byte segment: 255 minus the `2 * parity` parity bytes.
    #[inline]
    #[must_use]
    pub const fn segment_length(&self) -> u8 {
        255 - self.parity_bytes()
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};
    use crate::MAX_PARITY;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_segment_length_zero_parity() -> Result<(), TestError> {
        let header = LongEccHeader::new(0, OverlapFactor::Simple, 82, 50, 0)?;

        assert_eq!(header.segment_length(), 255);

        Ok(())
    }

    #[test]
    fn test_segment_length_subtracts_parity_bytes() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0)?;

        assert_eq!(header.segment_length(), 251);

        Ok(())
    }

    #[test]
    fn test_segment_length_maximum_parity() -> Result<(), TestError> {
        let header = LongEccHeader::new(MAX_PARITY, OverlapFactor::Simple, 208, 50, 0)?;

        assert_eq!(header.segment_length(), 129);

        Ok(())
    }
}
