use crate::LongEccDecodeError;

use super::checksums::xxh64;
use super::{LongEccHeader, HEADER_SIZE};

/// Fast validation using the codeword checksum.
///
/// Verifies the message and parity bytes without running error correction,
/// so corruption anywhere after the header is detected. A buffer whose length
/// differs from the header's full length is not a pristine codeword.
///
/// Returns the parsed header for a pristine codeword, `Ok(None)` for a
/// corrupted payload or a buffer too short to hold a header, and an error
/// if a full-sized header fails to parse.
///
/// # Errors
/// Returns an error if a full-sized header fails to parse, or if the header's
/// full length does not fit in `usize`.
pub fn fast_validate(codeword: &[u8]) -> Result<Option<LongEccHeader>, LongEccDecodeError> {
    if codeword.len() < HEADER_SIZE {
        return Ok(None);
    }

    let header = LongEccHeader::from_byte_slice(codeword)?;

    // A header that required correction means the codeword is not pristine.
    if codeword[..HEADER_SIZE] != header.to_bytes() {
        return Ok(None);
    }

    let full_length = usize::try_from(header.full_length())?;

    // Reject truncated codewords and codewords carrying trailing bytes.
    if codeword.len() != full_length {
        return Ok(None);
    }

    let Some(payload) = codeword.get(HEADER_SIZE..full_length) else {
        return Ok(None);
    };

    // Validate XXH64 of the message and parity
    Ok((xxh64(payload) == header.checksum()).then_some(header))
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use super::super::{encode, OverlapFactor};
    use super::fast_validate;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_fast_validate_pristine_codeword() -> Result<(), TestError> {
        let message = b"Pristine Codeword".to_buffer()?;

        let encoded = encode(&message, 2, OverlapFactor::Simple)?;

        assert!(fast_validate(&encoded)?.is_some());

        Ok(())
    }

    #[test]
    fn test_fast_validate_rejects_trailing_bytes() -> Result<(), TestError> {
        let message = b"Trailing Bytes".to_buffer()?;

        let encoded = encode(&message, 2, OverlapFactor::Simple)?;

        let mut extended = encoded.to_vec();

        extended.push(0xAB);

        assert!(fast_validate(&extended)?.is_none());

        Ok(())
    }
}
