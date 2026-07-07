use crate::long::LongEccHeader;

use super::super::OverlapFactor;

impl LongEccHeader {
    #[inline]
    #[must_use]
    pub const fn overlap_factor(&self) -> OverlapFactor {
        OverlapFactor::from_u8(self.parity)
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_overlap_factor_roundtrips_through_header() -> Result<(), TestError> {
        for factor in [
            OverlapFactor::Simple,
            OverlapFactor::Double,
            OverlapFactor::Triple,
            OverlapFactor::Quadruple,
        ] {
            let header = LongEccHeader::new(2, factor, 86, 50, 0)?;

            assert_eq!(header.overlap_factor(), factor);
        }

        Ok(())
    }
}
