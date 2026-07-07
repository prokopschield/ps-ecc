use crate::long::LongEccHeader;

impl LongEccHeader {
    #[inline]
    #[must_use]
    pub const fn header_parity(&self) -> [u8; 8] {
        self.header_parity
    }
}

#[cfg(test)]
mod tests {
    use crate::long::header::RS;
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_header_parity_matches_generated() -> Result<(), TestError> {
        let header = LongEccHeader::new(3, OverlapFactor::Triple, 88, 50, 7)?;

        let expected = RS.generate_parity(&header.to_bytes()[..24])?;

        assert_eq!(&header.header_parity()[..], &expected[..]);

        Ok(())
    }

    #[test]
    fn test_header_parity_roundtrips_through_to_bytes() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0)?;

        assert_eq!(header.to_bytes()[24..32], header.header_parity());

        Ok(())
    }
}
