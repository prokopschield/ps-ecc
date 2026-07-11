use ps_buffer::Buffer;

use crate::{Codeword, LongEccDecodeError};

use super::correct_in_place::correct_in_place;
use super::fast_validate;
use super::HEADER_SIZE;

/// Decodes a codeword to extract the original message.
pub fn decode(codeword: &[u8]) -> Result<Codeword<'_>, LongEccDecodeError> {
    let (header, is_ok) = fast_validate(codeword)?;

    if is_ok {
        return Ok(Codeword {
            codeword: codeword.into(),
            range: HEADER_SIZE..HEADER_SIZE + usize::try_from(header.message_length())?,
        });
    }

    let mut buffer = Buffer::from_slice(codeword)?;
    let header = correct_in_place(&mut buffer)?;

    let message_length = usize::try_from(header.message_length())?;
    let codeword = Codeword {
        codeword: buffer.into(),
        range: HEADER_SIZE..HEADER_SIZE + message_length,
    };

    Ok(codeword)
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::LongEccDecodeError;

    use super::super::{decode, encode, LongEccHeader, OverlapFactor, HEADER_SIZE};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_long_ecc_decode_no_errors() -> Result<(), TestError> {
        let message = b"Decode No Errors".to_buffer()?;
        let parity: u8 = 2;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_one_error() -> Result<(), TestError> {
        let message = b"Decode One Error".to_buffer()?;
        let parity: u8 = 2;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        encoded[HEADER_SIZE + 2] ^= 0b0000_0001;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_corrupted_parity() -> Result<(), TestError> {
        let message = b"Decode Corrupted Parity".to_buffer()?;
        let parity: u8 = 2;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        // Corrupt a parity byte; correction repairs it, so the checksum over
        // the message and parity verifies after correction.
        encoded[HEADER_SIZE + message.len()] ^= 0b0000_0001;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_too_many_errors() -> Result<(), TestError> {
        let message = b"Decode Too Many".to_buffer()?;
        let parity: u8 = 1;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        encoded[HEADER_SIZE + 1] ^= 0b0000_0001;
        encoded[HEADER_SIZE + 3] ^= 0b0000_0010;

        let result = decode(&encoded);

        assert!(matches!(
            result,
            Err(LongEccDecodeError::RSDecodeError(
                crate::RSDecodeError::RSComputeErrorsError(
                    crate::RSComputeErrorsError::TooManyErrors
                )
            ))
        ));

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_truncated_data() -> Result<(), TestError> {
        // Full length 44 covers the header (32), the message (10), and parity (2).
        let header = LongEccHeader::new(1, OverlapFactor::Simple, 44, 10, 0)?;

        let mut truncated = header.to_bytes().to_vec();

        truncated.extend_from_slice(&[0u8; 5]); // Only 5 message bytes instead of 10

        let result = decode(&truncated);

        assert!(matches!(result, Err(LongEccDecodeError::InvalidCodeword)));

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_discards_trailing_bytes() -> Result<(), TestError> {
        let message = b"Trailing Bytes".to_buffer()?;
        let parity: u8 = 2;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let mut extended = encoded.to_vec();

        extended.extend_from_slice(&[0xAB; 40]);

        let decoded = decode(&extended)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_incorrect_full_length() -> Result<(), TestError> {
        let message = b"Wrong Length".to_buffer()?;
        let parity: u8 = 1;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let mut incorrect_length = encoded.to_vec();

        // Corrupt a message byte so that correction is attempted, then truncate
        // the codeword so that the final segment's parity is incomplete.
        incorrect_length[HEADER_SIZE] ^= 0b0000_0001;
        incorrect_length.pop();

        let result = decode(&incorrect_length);

        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_zero_parity() -> Result<(), TestError> {
        let message = b"Zero Parity Decode".to_buffer()?;
        let parity: u8 = 0;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }
}
