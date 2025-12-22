use ps_buffer::Buffer;

use crate::{codeword::Codeword, long, DecodeError, EncodeError, ReedSolomon};

/// Encodes a message by adding an error-correcting code.
/// # Errors
/// - `RSConstructorError` is returned if `len(message) + 2 * parity` > `255`.
/// - `RSEncodeError` is returned if encoding fails for any reason.
pub fn encode(message: &[u8], parity: u8) -> Result<Buffer, EncodeError> {
    if message.len() + (usize::from(parity) << 1) > 0xff {
        let segment_length = 0xFF - (parity << 1);
        let codeword = long::encode(message, parity, segment_length, segment_length)?;

        return Ok(codeword);
    }

    let rs = ReedSolomon::new(parity)?;

    Ok(rs.encode(message)?)
}

/// Verifies the error-correcting code and returns the message.
/// # Errors
/// - `InputTooLarge` is returned if `len(received)` > 255 bytes.
/// - `InsufficientParityBytes` is returned if `parity > length / 2`.
/// - `RSDecodeError` is returned if decoding fails for any reason.
pub fn decode(received: &[u8], parity: u8) -> Result<Codeword<'_>, DecodeError> {
    if let Ok(length) = u8::try_from(received.len()) {
        if parity > length >> 1 {
            return Err(DecodeError::InsufficientParityBytes(parity, length));
        }

        let rs = ReedSolomon::new(parity)?;

        Ok(rs.decode(received)?)
    } else {
        Ok(long::decode(received)?)
    }
}

#[must_use]
/// Validates that a received codeword isn't corrupted.
pub fn validate(received: &[u8], parity: u8) -> bool {
    if let Ok(length) = u8::try_from(received.len()) {
        if parity > length >> 1 {
            return false;
        }

        let Ok(rs) = ReedSolomon::new(parity) else {
            return false;
        };

        match rs.validate(received) {
            Ok(None) => true,
            Ok(Some(_)) | Err(_) => false,
        }
    } else {
        long::fast_validate(received).unwrap_or_default()
    }
}

#[cfg(test)]
mod tests {
    use crate::EccError;

    use super::{decode, encode};

    #[test]
    fn ecc_works() -> Result<(), EccError> {
        let test_str = "Strč prst skrz krk! ¯\\_(ツ)_/¯".as_bytes();
        let mut encoded = encode(test_str, 13)?;

        for i in 0..13 {
            let index = (i * 37) % encoded.len();
            encoded[index] ^= (i * index + 13).to_le_bytes()[0];
            let decoded = decode(&encoded, 13)?;

            assert_eq!(test_str, &decoded[..]);
        }

        Ok(())
    }
}

#[cfg(test)]
mod validate_tests {
    use crate::{
        long, validate, LongEccConstructorError, LongEccDecodeError, LongEccEncodeError,
        LongEccToBytesError, RSEncodeError, ReedSolomon,
    };

    use ps_buffer::ToBuffer;

    #[derive(thiserror::Error, Debug)]
    enum TestError {
        #[error(transparent)]
        LongEccConstructor(#[from] LongEccConstructorError),
        #[error(transparent)]
        LongEccEncode(#[from] LongEccEncodeError),
        #[error(transparent)]
        LongEccDecode(#[from] LongEccDecodeError),
        #[error(transparent)]
        LongEccToBytes(#[from] LongEccToBytesError),
        #[error(transparent)]
        Buffer(#[from] ps_buffer::BufferError),
        #[error(transparent)]
        RSConstructorError(#[from] crate::RSConstructorError),
        #[error(transparent)]
        RSEncodeError(#[from] RSEncodeError),
    }

    #[test]
    fn test_validate_short_data_valid_no_errors() -> Result<(), TestError> {
        let data = b"test";
        let parity = 2;
        let rs = ReedSolomon::new(parity)?;
        let codeword = rs.encode(data)?;

        assert!(validate(&codeword, parity));
        Ok(())
    }

    #[test]
    fn test_validate_short_data_invalid_with_errors() -> Result<(), TestError> {
        let data = b"test";
        let parity = 2;
        let rs = ReedSolomon::new(parity)?;
        let mut codeword = rs.encode(data)?;

        // Introduce errors
        codeword[0] ^= 1;
        codeword[1] ^= 1;
        codeword[2] ^= 1;

        assert!(!validate(&codeword, parity));
        Ok(())
    }

    #[test]
    fn test_validate_short_data_parity_too_large() {
        let data = b"test"; // 4 bytes
        let parity = 3; // 3 > 4/2 (2) - parity is too large

        // Should return false when parity is too large
        assert!(!validate(data, parity));
    }

    #[test]
    fn test_validate_short_data_rs_constructor_error() {
        let data = b"test";
        let parity = 255; // Invalid parity value that should cause RS constructor to fail

        // Should return false when RS constructor fails
        assert!(!validate(data, parity));
    }

    #[test]
    fn test_validate_short_data_correctable_errors() -> Result<(), TestError> {
        let data = b"test";
        let parity = 2;
        let rs = ReedSolomon::new(parity)?;
        let mut codeword = rs.encode(data)?;

        // Introduce correctable errors (1 error with parity=2)
        codeword[0] ^= 1;

        // Validation should return false because there are errors (even if correctable)
        assert!(!validate(&codeword, parity));

        Ok(())
    }

    #[test]
    fn test_validate_long_data_valid_no_errors() -> Result<(), TestError> {
        let message = b"This is a longer message that will use long ECC".repeat(7);
        let parity = 2;
        let segment_length = 20;
        let segment_distance = 16;

        let encoded = long::encode(&message, parity, segment_length, segment_distance)?;

        // Valid data should pass validation
        assert!(validate(&encoded, parity));

        Ok(())
    }

    #[test]
    fn test_validate_long_data_invalid_with_errors() -> Result<(), TestError> {
        let message = b"This is a longer message that will use long ECC".to_buffer()?;
        let parity = 2;
        let segment_length = 20;
        let segment_distance = 16;

        let mut encoded = long::encode(&message, parity, segment_length, segment_distance)?;

        // Introduce errors in the message
        encoded[32] ^= 1;
        encoded[37] ^= 1;

        // Invalid data should fail validation
        assert!(!validate(&encoded, parity));
        Ok(())
    }

    #[test]
    fn test_validate_long_data_fast_path_valid() -> Result<(), TestError> {
        let message = b"Fast path validation test".repeat(12);
        let parity = 1;
        let segment_length = 15;
        let segment_distance = 12;

        let encoded = long::encode(&message, parity, segment_length, segment_distance)?;

        // Valid data should pass fast validation
        assert!(validate(&encoded, parity));

        Ok(())
    }

    #[test]
    fn test_validate_empty_data() {
        let data = b"";
        let parity = 0;

        // Empty data should be handled gracefully
        assert!(validate(data, parity));
    }

    #[test]
    fn test_validate_single_byte() -> Result<(), TestError> {
        let data = b"A";
        let parity = 1;
        let rs = ReedSolomon::new(parity)?;
        let codeword = rs.encode(data)?;

        assert!(validate(&codeword, parity));
        Ok(())
    }

    #[test]
    fn test_validate_large_short_data() -> Result<(), TestError> {
        let data = b"This is exactly 32 bytes of test data!!";
        let parity = 4;
        let rs = ReedSolomon::new(parity)?;
        let codeword = rs.encode(data)?;

        assert!(validate(&codeword, parity));
        Ok(())
    }

    #[test]
    fn test_validate_edge_case_parity_equals_length_div_2() -> Result<(), TestError> {
        let data = b"test"; // 4 bytes
        let parity = 2; // 2 == 4/2 - edge case

        let rs = ReedSolomon::new(parity)?;
        let codeword = rs.encode(data)?;

        assert!(validate(&codeword, parity));
        Ok(())
    }

    #[test]
    fn test_validate_edge_case_parity_just_over_length_div_2() {
        let data = b"test"; // 4 bytes
        let parity = 3; // 3 > 4/2 (2) - parity is too large

        // Should return false when parity is too large
        assert!(!validate(data, parity));
    }

    #[test]
    fn test_validate_long_data_with_zero_parity() -> Result<(), TestError> {
        let message = b"Zero parity test".to_buffer()?;
        let parity = 0;
        let segment_length = 15;
        let segment_distance = 12;

        let encoded = long::encode(&message, parity, segment_length, segment_distance)?;

        // Zero parity should still validate correctly
        assert!(validate(&encoded, parity));
        Ok(())
    }

    #[test]
    fn test_validate_long_data_header_corrupted() -> Result<(), TestError> {
        let message = b"Header corruption test".to_buffer()?;
        let parity = 2;
        let segment_length = 15;
        let segment_distance = 12;

        let mut encoded = long::encode(&message, parity, segment_length, segment_distance)?;

        // Corrupt header data
        encoded[0] ^= 1;
        encoded[5] ^= 1;

        // Corrupted header should fail validation
        assert!(!validate(&encoded, parity));
        Ok(())
    }

    #[test]
    fn test_validate_short_data_length_conversion_error() {
        let data: Vec<u8> = vec![0x42; 300]; // 300 bytes > 255
        let parity = 2;

        // This should fall back to long::fast_validate
        assert!(!validate(&data, parity));
    }

    #[test]
    fn test_validate_short_data_unrecoverable_errors() -> Result<(), TestError> {
        let data = b"test data";
        let parity = 1;
        let rs = ReedSolomon::new(parity)?;
        let mut codeword = rs.encode(data)?;

        // Introduce unrecoverable errors (more errors than parity can correct)
        codeword[0] ^= 1;
        codeword[1] ^= 1; // 2 errors with parity=1 - unrecoverable

        assert!(!validate(&codeword, parity));
        Ok(())
    }

    #[test]
    fn test_validate_long_data_corrupted_parity() -> Result<(), TestError> {
        let message = b"Corrupted parity test".to_buffer()?;
        let parity = 2;
        let segment_length = 15;
        let segment_distance = 12;

        let mut encoded = long::encode(&message, parity, segment_length, segment_distance)?;

        // Corrupt parity bytes
        let parity_start = 32 + message.len();
        if parity_start < encoded.len() {
            encoded[parity_start] ^= 1;
        }

        // Corrupted parity should fail validation
        assert!(!validate(&encoded, parity));
        Ok(())
    }
}
