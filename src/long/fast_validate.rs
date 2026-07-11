use crate::LongEccFastValidateError;

use super::checksums::xxh64;
use super::{LongEccHeader, HEADER_SIZE};

/// Fast validation using the codeword checksum.
///
/// Verifies the message and parity bytes without running error correction,
/// so corruption anywhere after the header is detected. A buffer whose length
/// differs from the header's full length is not a pristine codeword.
///
/// Returns the parsed header together with a flag that is `true` only for a
/// pristine codeword: one whose stored header bytes needed no correction,
/// whose length matches the header's full length, and whose payload matches
/// the header's checksum.
///
/// # Errors
/// Returns an error if the header fails to parse, including when the buffer
/// is too short to hold a header, or if the header's full length does not
/// fit in `usize`.
pub fn fast_validate(codeword: &[u8]) -> Result<(LongEccHeader, bool), LongEccFastValidateError> {
    let header = LongEccHeader::from_byte_slice(codeword)?;

    // A header that required correction means the codeword is not pristine.
    if codeword[..HEADER_SIZE] != header.to_bytes() {
        return Ok((header, false));
    }

    let full_length = usize::try_from(header.full_length())?;

    // Reject truncated codewords and codewords carrying trailing bytes.
    if codeword.len() != full_length {
        return Ok((header, false));
    }

    let Some(payload) = codeword.get(HEADER_SIZE..full_length) else {
        return Ok((header, false));
    };

    // Validate XXH64 of the message and parity
    Ok((header, xxh64(payload) == header.checksum()))
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

        assert!(fast_validate(&encoded)?.1);

        Ok(())
    }

    #[test]
    fn test_fast_validate_rejects_trailing_bytes() -> Result<(), TestError> {
        let message = b"Trailing Bytes".to_buffer()?;

        let encoded = encode(&message, 2, OverlapFactor::Simple)?;

        let mut extended = encoded.to_vec();

        extended.push(0xAB);

        assert!(!fast_validate(&extended)?.1);

        Ok(())
    }
}
