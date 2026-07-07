use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Stride between consecutive segment starts: `segment_length` divided by the overlap
    /// factor, so segments cover each byte `overlap_factor` times.
    #[inline]
    #[must_use]
    pub const fn segment_distance(&self) -> u8 {
        self.segment_length() / self.overlap_factor().count()
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_segment_distance_divides_by_overlap_count() -> Result<(), TestError> {
        // With parity 2, the segment length is 251; the stride is 251 divided
        // by the overlap count, rounded down.
        let cases = [
            (OverlapFactor::Simple, 251),
            (OverlapFactor::Double, 125),
            (OverlapFactor::Triple, 83),
            (OverlapFactor::Quadruple, 62),
        ];

        for (factor, expected) in cases {
            let header = LongEccHeader::new(2, factor, 86, 50, 0)?;

            assert_eq!(header.segment_distance(), expected);
        }

        Ok(())
    }

    #[test]
    fn test_segment_distance_equals_segment_length_without_overlap() -> Result<(), TestError> {
        let header = LongEccHeader::new(4, OverlapFactor::Simple, 90, 50, 0)?;

        assert_eq!(header.segment_distance(), header.segment_length());

        Ok(())
    }
}
