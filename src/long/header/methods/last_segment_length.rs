use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Data length of the final ECC segment, always shorter than `segment_length`.
    ///
    /// The constructors enforce `2 * parity < segment_distance`, so the
    /// subtraction cannot underflow.
    #[inline]
    #[must_use]
    // The remainder is always shorter than `segment_length`, so it fits in `u8`.
    #[allow(clippy::cast_possible_truncation)]
    pub const fn last_segment_length(&self) -> u8 {
        let new_bytes_per_segment = (self.segment_distance() - self.parity_bytes()) as u32;
        let full_segments = self.segment_count() - 1;
        let covered = full_segments.saturating_mul(new_bytes_per_segment);

        self.message_length().saturating_sub(covered) as u8
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_last_segment_length_single_segment() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 136, 100, 0)?;

        assert_eq!(header.last_segment_length(), 100);

        Ok(())
    }

    #[test]
    fn test_last_segment_length_empty_message() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 36, 0, 0)?;

        assert_eq!(header.last_segment_length(), 0);

        Ok(())
    }

    #[test]
    fn test_last_segment_length_multiple_segments() -> Result<(), TestError> {
        // Three segments cover 600 bytes; the first two contribute 239 new
        // bytes each, leaving 600 - 2 * 239 = 122 bytes for the last segment.
        let header = LongEccHeader::new(4, OverlapFactor::Simple, 656, 600, 0)?;

        assert_eq!(header.last_segment_length(), 122);
        assert!(header.last_segment_length() < header.segment_length());

        Ok(())
    }

    #[test]
    fn test_last_segment_length_boundary_spill() -> Result<(), TestError> {
        // A 247-byte message with parity 4 spills exactly 247 - 239 = 8 bytes
        // into the second segment.
        let header = LongEccHeader::new(4, OverlapFactor::Simple, 295, 247, 0)?;

        assert_eq!(header.last_segment_length(), 8);

        Ok(())
    }
}
