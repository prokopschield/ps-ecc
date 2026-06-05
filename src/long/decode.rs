use ps_buffer::Buffer;

use crate::{Codeword, LongEccDecodeError};

use super::correct_in_place::correct_in_place;
use super::HEADER_SIZE;

/// Decode a codeword to extract the original message
pub fn decode(codeword: &[u8]) -> Result<Codeword<'_>, LongEccDecodeError> {
    let mut buffer = Buffer::from_slice(codeword)?;
    let header = correct_in_place(&mut buffer)?;

    let message_length = usize::try_from(header.message_length)?;
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

    use super::super::{decode, encode, LongEccHeader, HEADER_SIZE};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_long_ecc_decode_no_errors() -> Result<(), TestError> {
        let message = b"Decode No Errors".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_one_error() -> Result<(), TestError> {
        let message = b"Decode One Error".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;

        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;

        encoded[HEADER_SIZE + 2] ^= 0b0000_0001;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_too_many_errors() -> Result<(), TestError> {
        let message = b"Decode Too Many".to_buffer()?;
        let parity: u8 = 1;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;

        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;

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
        let header = LongEccHeader {
            full_length: 20,
            message_length: 10,
            parity: 1,
            segment_length: 8,
            segment_distance: 6,
            last_segment_length: 8,
            crc32: 0,
            xxh64: 0,
        };

        let mut truncated = header.to_bytes()?.to_vec();

        truncated.extend_from_slice(&[0u8; 5]); // Only 5 message bytes instead of 10

        let result = decode(&truncated);

        eprintln!("{result:?}");

        assert!(matches!(
            result,
            Err(crate::LongEccDecodeError::RSDecodeError(
                crate::RSDecodeError::RSComputeErrorsError(
                    crate::RSComputeErrorsError::TooManyErrors
                )
            ))
        ));

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_incorrect_full_length() -> Result<(), TestError> {
        let message = b"Wrong Length".to_buffer()?;
        let parity: u8 = 1;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let mut incorrect_length = encoded.to_vec();

        incorrect_length.pop();

        let result = decode(&incorrect_length);

        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_zero_parity() -> Result<(), TestError> {
        let message = b"Zero Parity Decode".to_buffer()?;
        let parity: u8 = 0;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }
}
