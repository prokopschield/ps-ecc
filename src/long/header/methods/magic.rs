use crate::long::LongEccHeader;

impl LongEccHeader {
    #[inline]
    #[must_use]
    pub const fn magic(&self) -> u16 {
        self.magic
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor, LONG_ECC_HEADER_MAGIC};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_magic_returns_header_magic() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0)?;

        assert_eq!(header.magic(), LONG_ECC_HEADER_MAGIC);

        Ok(())
    }

    #[test]
    fn test_magic_is_ascii_lh() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0)?;

        assert_eq!(header.magic().to_be_bytes(), *b"LH");

        Ok(())
    }
}
