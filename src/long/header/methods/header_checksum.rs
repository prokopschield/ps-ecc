use crate::long::LongEccHeader;

impl LongEccHeader {
    #[inline]
    #[must_use]
    pub const fn header_checksum(&self) -> u32 {
        self.header_checksum
    }
}

#[cfg(test)]
mod tests {
    use crate::long::header::utils::calculate_header_checksum;
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_header_checksum_matches_calculated() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Double, 86, 50, 0)?;

        let packed_parity = 2 | OverlapFactor::Double.to_u8();
        let expected = calculate_header_checksum(header.version(), packed_parity, 86, 50);

        assert_eq!(header.header_checksum(), expected);

        Ok(())
    }

    #[test]
    fn test_header_checksum_differs_across_headers() -> Result<(), TestError> {
        let first = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0)?;
        let second = LongEccHeader::new(2, OverlapFactor::Simple, 87, 51, 0)?;

        assert_ne!(first.header_checksum(), second.header_checksum());

        Ok(())
    }
}
