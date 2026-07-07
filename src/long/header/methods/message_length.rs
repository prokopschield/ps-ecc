use crate::long::LongEccHeader;

impl LongEccHeader {
    #[inline]
    #[must_use]
    pub const fn message_length(&self) -> u32 {
        self.message_length
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_message_length_returns_constructor_value() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0)?;

        assert_eq!(header.message_length(), 50);

        Ok(())
    }

    #[test]
    fn test_message_length_zero() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 36, 0, 0)?;

        assert_eq!(header.message_length(), 0);

        Ok(())
    }
}
