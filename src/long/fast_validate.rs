use crate::LongEccDecodeError;

use super::checksums::{crc32, xxh64};
use super::{LongEccHeader, HEADER_SIZE};

/// Fast validation using checksums.
///
/// Returns the parsed header if the codeword passes all checksum checks,
/// and `Ok(None)` if any check fails.
///
/// # Errors
/// Returns an error if the header cannot be parsed.
pub fn fast_validate(codeword: &[u8]) -> Result<Option<LongEccHeader>, LongEccDecodeError> {
    if codeword.len() < HEADER_SIZE {
        return Ok(None);
    }

    let header = LongEccHeader::from_bytes(codeword)?;

    // Validate CRC32 of message
    let message_start = HEADER_SIZE;
    let message_end = HEADER_SIZE + header.message_length as usize;

    if message_end > codeword.len() {
        return Ok(None);
    }

    let message = &codeword[message_start..message_end];
    let calculated_crc32 = crc32(message);

    if calculated_crc32 != header.crc32 {
        return Ok(None);
    }

    // Validate XXH64 of message + parity

    let calculated_xxh64 = xxh64(&codeword[HEADER_SIZE..]);

    if calculated_xxh64 != header.xxh64 {
        return Ok(None);
    }

    Ok(Some(header))
}
