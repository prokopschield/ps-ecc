use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Returns the XXH64 checksum of the message and parity bytes.
    #[inline]
    #[must_use]
    pub const fn checksum(&self) -> u64 {
        self.checksum
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_checksum_returns_constructor_value() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0xDEAD_BEEF_CAFE_BABE)?;

        assert_eq!(header.checksum(), 0xDEAD_BEEF_CAFE_BABE);

        Ok(())
    }

    #[test]
    fn test_checksum_maximum_value() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, u64::MAX)?;

        assert_eq!(header.checksum(), u64::MAX);

        Ok(())
    }
}
