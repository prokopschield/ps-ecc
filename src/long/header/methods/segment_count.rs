use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Number of segments: every segment before the last contributes
    /// `segment_distance - 2 * parity` new data bytes, and the final, shorter
    /// segment covers the remainder.
    ///
    /// The constructors enforce `2 * parity < segment_distance`, so the
    /// subtraction cannot underflow.
    #[inline]
    #[must_use]
    pub const fn segment_count(&self) -> u32 {
        let new_bytes_per_segment = (self.segment_distance() - self.parity_bytes()) as u32;
        let remaining_bytes = self
            .message_length()
            .saturating_sub(self.segment_length() as u32 - 1);

        remaining_bytes
            .div_ceil(new_bytes_per_segment)
            .saturating_add(1)
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_segment_count_single_segment() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 136, 100, 0)?;

        assert_eq!(header.segment_count(), 1);

        Ok(())
    }

    #[test]
    fn test_segment_count_empty_message() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 36, 0, 0)?;

        assert_eq!(header.segment_count(), 1);

        Ok(())
    }

    #[test]
    fn test_segment_count_multiple_segments() -> Result<(), TestError> {
        // With parity 4, a segment holds 247 data bytes and each full segment
        // contributes 239 new bytes, so 600 bytes require three segments.
        let header = LongEccHeader::new(4, OverlapFactor::Simple, 656, 600, 0)?;

        assert_eq!(header.segment_count(), 3);

        Ok(())
    }

    #[test]
    fn test_segment_count_boundary_at_segment_length() -> Result<(), TestError> {
        // With parity 4, a segment holds 247 data bytes: a 246-byte message
        // fits in one segment, and a 247-byte message spills into a second.
        let single = LongEccHeader::new(4, OverlapFactor::Simple, 286, 246, 0)?;
        let double = LongEccHeader::new(4, OverlapFactor::Simple, 295, 247, 0)?;

        assert_eq!(single.segment_count(), 1);
        assert_eq!(double.segment_count(), 2);

        Ok(())
    }
}
