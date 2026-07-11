use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Returns the length of the full codeword, including the header and parity.
    #[inline]
    #[must_use]
    pub const fn full_length(&self) -> u32 {
        self.full_length
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_full_length_returns_constructor_value() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0)?;

        assert_eq!(header.full_length(), 86);

        Ok(())
    }

    #[test]
    fn test_full_length_maximum_value() -> Result<(), TestError> {
        // With zero parity, the codeword is the header plus the message, so the
        // maximum full length pairs with a message of `u32::MAX - 32` bytes.
        let header = LongEccHeader::new(0, OverlapFactor::Simple, u32::MAX, u32::MAX - 32, 0)?;

        assert_eq!(header.full_length(), u32::MAX);

        Ok(())
    }
}
