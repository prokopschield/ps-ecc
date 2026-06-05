use crate::LongEccDecodeError;

use super::checksums::{crc32, xxh64};
use super::{LongEccHeader, HEADER_SIZE};

/// Fast validation using checksums
pub fn fast_validate(codeword: &[u8]) -> Result<bool, LongEccDecodeError> {
    if codeword.len() < HEADER_SIZE {
        return Ok(false);
    }

    let header = LongEccHeader::from_bytes(codeword)?;

    // Validate CRC32 of message
    let message_start = HEADER_SIZE;
    let message_end = HEADER_SIZE + header.message_length as usize;

    if message_end > codeword.len() {
        return Ok(false);
    }

    let message = &codeword[message_start..message_end];
    let calculated_crc32 = crc32(message);

    if calculated_crc32 != header.crc32 {
        return Ok(false);
    }

    // Validate XXH64 of message + parity

    let calculated_xxh64 = xxh64(&codeword[HEADER_SIZE..]);

    if calculated_xxh64 != header.xxh64 {
        return Ok(false);
    }

    Ok(true)
}
