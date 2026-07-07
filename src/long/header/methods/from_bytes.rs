use ps_util::Array;

use crate::{
    long::{
        header::{
            magic::{LONG_ECC_HEADER_MAGIC, LONG_ECC_HEADER_VERSION},
            utils::calculate_header_checksum,
        },
        LongEccHeader, HEADER_SIZE,
    },
    RSDecodeError, ReedSolomon,
};

#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccHeaderFromBytesError {
    #[error("Incorrect magic number: {0:x}")]
    IncorrectMagic(u16),

    #[error("Incorrect version number: {0}")]
    InvalidVersion(u8),

    #[error("Header checksum incorrect.")]
    IncorrectChecksum,

    #[error("Header error correction failed: {0}")]
    CorrectionFailed(#[from] RSDecodeError),

    #[error("Message length {0} does not fit within full length {1}")]
    InvalidMessageLength(u32, u32),

    #[error("Invalid segment-to-parity ratio: {0} <= 2 * {1}")]
    InvalidSegmentParityRatio(u8, u8),

    #[error("Full length {0} does not match the derived codeword length {1}")]
    InvalidFullLength(u32, u64),
}

impl LongEccHeader {
    /// Parses a header from its serialized form, correcting errors via the header parity.
    ///
    /// # Errors
    /// Returns an error if header correction fails, or if the corrected header
    /// carries an incorrect magic number, version, or checksum, a message
    /// length that does not fit within the full length, an invalid
    /// segment-to-parity ratio, or a full length that does not match the
    /// derived codeword length.
    pub fn from_bytes(mut bytes: [u8; 32]) -> Result<Self, LongEccHeaderFromBytesError> {
        let (body, parity) = bytes.split_at_mut(24);

        ReedSolomon::correct_detached_in_place(parity, body)?;

        let magic = u16::from_be_bytes(*bytes.subarray(0));
        let version = bytes[2];
        let parity = bytes[3];
        let full_length = u32::from_be_bytes(*bytes.subarray(4));
        let message_length = u32::from_be_bytes(*bytes.subarray(8));
        let header_checksum = u32::from_be_bytes(*bytes.subarray(12));
        let checksum = u64::from_be_bytes(*bytes.subarray(16));
        let header_parity = *bytes.subarray(24);

        if magic != LONG_ECC_HEADER_MAGIC {
            return Err(LongEccHeaderFromBytesError::IncorrectMagic(magic));
        }

        if version != LONG_ECC_HEADER_VERSION {
            return Err(LongEccHeaderFromBytesError::InvalidVersion(version));
        }

        if header_checksum != {
            calculate_header_checksum(version, parity, full_length, message_length)
        } {
            return Err(LongEccHeaderFromBytesError::IncorrectChecksum);
        }

        // The codeword must hold the header and the message; parity only adds
        // to it. The exact-length check below subsumes this bound, but this
        // check runs first and is kept for its more specific
        // `InvalidMessageLength` diagnostic.
        if u64::from(full_length) < u64::from(message_length) + HEADER_SIZE as u64 {
            return Err(LongEccHeaderFromBytesError::InvalidMessageLength(
                message_length,
                full_length,
            ));
        }

        let header = Self {
            magic,
            version,
            parity,
            full_length,
            message_length,
            header_checksum,
            checksum,
            header_parity,
        };

        // Reject headers whose parity leaves no room for new data within a segment;
        // the derived segment methods assume `2 * parity < segment_distance`.
        if header.parity_bytes() >= header.segment_distance() {
            return Err(LongEccHeaderFromBytesError::InvalidSegmentParityRatio(
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
            return Err(LongEccHeaderFromBytesError::InvalidFullLength(
                full_length,
                expected_length,
            ));
        }

        Ok(header)
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{
        header::{
            magic::{LONG_ECC_HEADER_MAGIC, LONG_ECC_HEADER_VERSION},
            utils::calculate_header_checksum,
            RS,
        },
        LongEccHeader, OverlapFactor,
    };

    use super::LongEccHeaderFromBytesError;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_from_bytes_rejects_invalid_segment_parity_ratio() -> Result<(), TestError> {
        // Craft a header that is valid on the wire (magic, version, checksum, and
        // RS parity all check out) but carries an invalid geometry: with parity 63
        // and quadruple overlap, each segment holds 126 parity bytes, but
        // consecutive segments start only (255 - 126) / 4 = 32 bytes apart.
        let parity = 63 | OverlapFactor::Quadruple.to_u8();
        let full_length = 1000u32;
        let message_length = 500u32;
        let header_checksum =
            calculate_header_checksum(LONG_ECC_HEADER_VERSION, parity, full_length, message_length);

        let mut bytes = [0u8; 32];

        bytes[..2].copy_from_slice(&LONG_ECC_HEADER_MAGIC.to_be_bytes());
        bytes[2] = LONG_ECC_HEADER_VERSION;
        bytes[3] = parity;
        bytes[4..8].copy_from_slice(&full_length.to_be_bytes());
        bytes[8..12].copy_from_slice(&message_length.to_be_bytes());
        bytes[12..16].copy_from_slice(&header_checksum.to_be_bytes());

        let header_parity = RS.generate_parity(&bytes[..24])?;

        bytes[24..].copy_from_slice(&header_parity);

        let result = LongEccHeader::from_bytes(bytes);

        assert_eq!(
            result,
            Err(LongEccHeaderFromBytesError::InvalidSegmentParityRatio(
                32, 63
            ))
        );

        Ok(())
    }

    #[test]
    fn test_from_bytes_rejects_message_length_exceeding_full_length() -> Result<(), TestError> {
        // Craft a header that is valid on the wire but whose message cannot fit:
        // full length 100 leaves room for only 68 message bytes after the header.
        let parity = 2 | OverlapFactor::Simple.to_u8();
        let full_length = 100u32;
        let message_length = 100u32;
        let header_checksum =
            calculate_header_checksum(LONG_ECC_HEADER_VERSION, parity, full_length, message_length);

        let mut bytes = [0u8; 32];

        bytes[..2].copy_from_slice(&LONG_ECC_HEADER_MAGIC.to_be_bytes());
        bytes[2] = LONG_ECC_HEADER_VERSION;
        bytes[3] = parity;
        bytes[4..8].copy_from_slice(&full_length.to_be_bytes());
        bytes[8..12].copy_from_slice(&message_length.to_be_bytes());
        bytes[12..16].copy_from_slice(&header_checksum.to_be_bytes());

        let header_parity = RS.generate_parity(&bytes[..24])?;

        bytes[24..].copy_from_slice(&header_parity);

        let result = LongEccHeader::from_bytes(bytes);

        assert_eq!(
            result,
            Err(LongEccHeaderFromBytesError::InvalidMessageLength(100, 100))
        );

        Ok(())
    }

    #[test]
    fn test_from_bytes_rejects_incorrect_full_length() -> Result<(), TestError> {
        // Craft a header that is valid on the wire but inflated: with parity 1,
        // the codeword spans the 32-byte header, 10 message bytes, and 2 parity
        // bytes, so the full length must be 44, not 100.
        let parity = 1 | OverlapFactor::Simple.to_u8();
        let full_length = 100u32;
        let message_length = 10u32;
        let header_checksum =
            calculate_header_checksum(LONG_ECC_HEADER_VERSION, parity, full_length, message_length);

        let mut bytes = [0u8; 32];

        bytes[..2].copy_from_slice(&LONG_ECC_HEADER_MAGIC.to_be_bytes());
        bytes[2] = LONG_ECC_HEADER_VERSION;
        bytes[3] = parity;
        bytes[4..8].copy_from_slice(&full_length.to_be_bytes());
        bytes[8..12].copy_from_slice(&message_length.to_be_bytes());
        bytes[12..16].copy_from_slice(&header_checksum.to_be_bytes());

        let header_parity = RS.generate_parity(&bytes[..24])?;

        bytes[24..].copy_from_slice(&header_parity);

        let result = LongEccHeader::from_bytes(bytes);

        assert_eq!(
            result,
            Err(LongEccHeaderFromBytesError::InvalidFullLength(100, 44))
        );

        Ok(())
    }

    #[test]
    fn test_from_bytes_rejects_zero_checksum_for_zero_parity() -> Result<(), TestError> {
        // The former multiplicative checksum was 0 for every zero-parity header,
        // validating arbitrary length fields; a stored 0 must now be rejected.
        let parity = OverlapFactor::Simple.to_u8();
        let full_length = 100u32;
        let message_length = 50u32;

        let mut bytes = [0u8; 32];

        bytes[..2].copy_from_slice(&LONG_ECC_HEADER_MAGIC.to_be_bytes());
        bytes[2] = LONG_ECC_HEADER_VERSION;
        bytes[3] = parity;
        bytes[4..8].copy_from_slice(&full_length.to_be_bytes());
        bytes[8..12].copy_from_slice(&message_length.to_be_bytes());

        let header_parity = RS.generate_parity(&bytes[..24])?;

        bytes[24..].copy_from_slice(&header_parity);

        let result = LongEccHeader::from_bytes(bytes);

        assert_eq!(result, Err(LongEccHeaderFromBytesError::IncorrectChecksum));

        Ok(())
    }
}
