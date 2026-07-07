use crate::long::LongEccHeader;

impl LongEccHeader {
    #[inline]
    #[must_use]
    pub const fn version(&self) -> u8 {
        self.version
    }
}

#[cfg(test)]
mod tests {
    use crate::long::header::magic::LONG_ECC_HEADER_VERSION;
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_version_is_current() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 0)?;

        assert_eq!(header.version(), LONG_ECC_HEADER_VERSION);
        assert_eq!(header.version(), 1);

        Ok(())
    }
}
