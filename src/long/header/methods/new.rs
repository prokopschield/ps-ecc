use crate::{
    long::{
        header::{
            magic::LONG_ECC_HEADER_VERSION, utils::calculate_header_checksum, LONG_ECC_HEADER_MAGIC,
        },
        LongEccHeader, OverlapFactor, HEADER_SIZE,
    },
    RSGenerateParityError, MAX_PARITY,
};

use super::super::RS;

/// Errors returned by [`LongEccHeader::new`](crate::LongEccHeader::new).
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccHeaderConstructorError {
    /// Propagated from generating the header parity.
    #[error("Generating parity failed: {0}")]
    GenerateParity(#[from] RSGenerateParityError),
    /// The full length does not match the codeword length derived from the
    /// message length and segment geometry.
    #[error("Full length {0} does not match the derived codeword length {1}.")]
    InvalidFullLength(u32, u64),
    /// The header and message do not fit within the full length.
    #[error("Message length {0} does not fit within full length {1}.")]
    InvalidMessageLength(u32, u32),
    /// The error-correction capability exceeds
    /// [`MAX_PARITY`](crate::MAX_PARITY).
    #[error("Invalid parity count: {0}.")]
    InvalidParityCount(u8),
    /// The parity bytes leave no room for new data within a segment.
    #[error("Invalid segment-to-parity ratio: {0} <= 2 * {1}.")]
    InvalidSegmentParityRatio(u8, u8),
}

impl LongEccHeader {
    /// Constructs a header, computing its checksum and parity.
    ///
    /// # Errors
    /// Returns `InvalidParityCount` if `parity` exceeds [`MAX_PARITY`],
    /// `InvalidMessageLength` if the header and message do not fit within
    /// `full_length`, `InvalidSegmentParityRatio` if the parity bytes leave
    /// no room for new data within a segment, `InvalidFullLength` if
    /// `full_length` does not match the codeword length derived from the
    /// message length and segment geometry, or `GenerateParity` if
    /// generating the header parity fails.
    pub fn new(
        parity: u8,
        overlap_factor: OverlapFactor,
        full_length: u32,
        message_length: u32,
        checksum: u64,
    ) -> Result<Self, LongEccHeaderConstructorError> {
        let magic = LONG_ECC_HEADER_MAGIC;
        let version = LONG_ECC_HEADER_VERSION;

        if parity > MAX_PARITY {
            return Err(LongEccHeaderConstructorError::InvalidParityCount(parity));
        }

        // The codeword must hold the header and the message; parity only adds
        // to it. The exact-length check below subsumes this bound, but this
        // check runs first and is kept for its more specific
        // `InvalidMessageLength` diagnostic.
        if u64::from(full_length) < u64::from(message_length) + HEADER_SIZE as u64 {
            return Err(LongEccHeaderConstructorError::InvalidMessageLength(
                message_length,
                full_length,
            ));
        }

        // Pack the overlap factor into the high two bits of the parity byte.
        let parity = parity | overlap_factor.to_u8();

        let mut header = Self {
            magic,
            version,
            parity,
            full_length,
            message_length,
            header_checksum: calculate_header_checksum(
                version,
                parity,
                full_length,
                message_length,
            ),
            checksum,
            header_parity: [0; 8],
        };

        // Reject headers whose parity leaves no room for new data within a segment;
        // the derived segment methods assume `2 * parity < segment_distance`.
        if header.parity_bytes() >= header.segment_distance() {
            return Err(LongEccHeaderConstructorError::InvalidSegmentParityRatio(
                header.segment_distance(),
                header.parity(),
            ));
        }

        // Defensively require the full length to match the codeword geometry
        // exactly: the header, the message, and the parity of every segment.
        let expected_length = HEADER_SIZE as u64
            + u64::from(message_length)
            + u64::from(header.parity_bytes()) * u64::from(header.segment_count());

        if u64::from(full_length) != expected_length {
            return Err(LongEccHeaderConstructorError::InvalidFullLength(
                full_length,
                expected_length,
            ));
        }

        let parity = RS.generate_parity(&header.to_bytes()[..24])?;

        header.header_parity.copy_from_slice(&parity);

        Ok(header)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        long::{LongEccHeader, OverlapFactor},
        MAX_PARITY,
    };

    use super::LongEccHeaderConstructorError;

    #[test]
    fn test_new_rejects_invalid_segment_parity_ratio() {
        // With parity 63 and quadruple overlap, each segment holds 126 parity bytes,
        // but consecutive segments start only (255 - 126) / 4 = 32 bytes apart.
        let result = LongEccHeader::new(63, OverlapFactor::Quadruple, 100, 50, 0);

        assert_eq!(
            result,
            Err(LongEccHeaderConstructorError::InvalidSegmentParityRatio(
                32, 63
            ))
        );
    }

    #[test]
    fn test_new_quadruple_overlap_boundary() {
        // With parity 25, the 50 parity bytes fit below the segment distance
        // of (255 - 50) / 4 = 51; with parity 26, the 52 parity bytes exceed
        // the segment distance of (255 - 52) / 4 = 50.
        assert!(LongEccHeader::new(25, OverlapFactor::Quadruple, 132, 50, 0).is_ok());
        assert_eq!(
            LongEccHeader::new(26, OverlapFactor::Quadruple, 132, 50, 0),
            Err(LongEccHeaderConstructorError::InvalidSegmentParityRatio(
                50, 26
            ))
        );
    }

    #[test]
    fn test_new_double_overlap_boundary() {
        // With parity 42, the 84 parity bytes fit below the segment distance
        // of (255 - 84) / 2 = 85; with parity 43, the 86 parity bytes exceed
        // the segment distance of (255 - 86) / 2 = 84.
        assert!(LongEccHeader::new(42, OverlapFactor::Double, 166, 50, 0).is_ok());
        assert_eq!(
            LongEccHeader::new(43, OverlapFactor::Double, 166, 50, 0),
            Err(LongEccHeaderConstructorError::InvalidSegmentParityRatio(
                84, 43
            ))
        );
    }

    #[test]
    fn test_new_triple_overlap_boundary() {
        // With parity 31, the 62 parity bytes fit below the segment distance
        // of (255 - 62) / 3 = 64; with parity 32, the 64 parity bytes exceed
        // the segment distance of (255 - 64) / 3 = 63.
        assert!(LongEccHeader::new(31, OverlapFactor::Triple, 144, 50, 0).is_ok());
        assert_eq!(
            LongEccHeader::new(32, OverlapFactor::Triple, 144, 50, 0),
            Err(LongEccHeaderConstructorError::InvalidSegmentParityRatio(
                63, 32
            ))
        );
    }

    #[test]
    fn test_new_rejects_message_length_exceeding_full_length() {
        // Full length 41 is one byte short of the 32-byte header plus 10 message bytes.
        let result = LongEccHeader::new(1, OverlapFactor::Simple, 41, 10, 0);

        assert_eq!(
            result,
            Err(LongEccHeaderConstructorError::InvalidMessageLength(10, 41))
        );
    }

    #[test]
    fn test_new_rejects_incorrect_full_length() {
        // With parity 1, the codeword spans the 32-byte header, 10 message
        // bytes, and 2 parity bytes, so the full length must be exactly 44.
        assert!(LongEccHeader::new(1, OverlapFactor::Simple, 44, 10, 0).is_ok());
        assert_eq!(
            LongEccHeader::new(1, OverlapFactor::Simple, 45, 10, 0),
            Err(LongEccHeaderConstructorError::InvalidFullLength(45, 44))
        );
    }

    #[test]
    fn test_new_accepts_max_parity_without_overlap() {
        // With maximum parity and no overlap, the 126 parity bytes fit below
        // the segment distance of 255 - 126 = 129.
        assert!(LongEccHeader::new(MAX_PARITY, OverlapFactor::Simple, 208, 50, 0).is_ok());
    }
}
