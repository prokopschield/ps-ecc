use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Returns the error-correction capability per segment: the number of byte errors each segment can recover from.
    #[inline]
    #[must_use]
    pub const fn parity(&self) -> u8 {
        self.parity & 0x3F
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};
    use crate::MAX_PARITY;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_parity_masks_overlap_factor() -> Result<(), TestError> {
        for factor in [
            OverlapFactor::Simple,
            OverlapFactor::Double,
            OverlapFactor::Triple,
            OverlapFactor::Quadruple,
        ] {
            let header = LongEccHeader::new(5, factor, 92, 50, 0)?;

            assert_eq!(header.parity(), 5);
        }

        Ok(())
    }

    #[test]
    fn test_parity_zero() -> Result<(), TestError> {
        let header = LongEccHeader::new(0, OverlapFactor::Simple, 82, 50, 0)?;

        assert_eq!(header.parity(), 0);

        Ok(())
    }

    #[test]
    fn test_parity_maximum() -> Result<(), TestError> {
        let header = LongEccHeader::new(MAX_PARITY, OverlapFactor::Simple, 208, 50, 0)?;

        assert_eq!(header.parity(), MAX_PARITY);

        Ok(())
    }
}
