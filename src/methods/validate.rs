use crate::{long, ReedSolomon};

/// Validates that a received codeword is pristine.
///
/// Returns `true` only if the codeword is entirely uncorrupted; for codewords
/// longer than 255 bytes, this includes carrying no bytes beyond the full
/// length recorded in the header. A `true` result implies that [`decode`](crate::decode)
/// succeeds. The converse does not hold: [`decode`](crate::decode) repairs correctable
/// corruption and discards trailing bytes, so it accepts input that this
/// function rejects.
#[must_use]
pub fn validate(received: &[u8], parity: u8) -> bool {
    if let Ok(length) = u8::try_from(received.len()) {
        if parity > length >> 1 {
            return false;
        }

        let Ok(rs) = ReedSolomon::new(parity) else {
            return false;
        };

        rs.validate(received).is_none()
    } else {
        matches!(long::fast_validate(received), Ok(Some(_)))
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        long, validate, LongEccDecodeError, LongEccEncodeError, RSEncodeError, ReedSolomon,
    };

    use ps_buffer::ToBuffer;

    #[derive(thiserror::Error, Debug)]
    enum TestError {
        #[error(transparent)]
        LongEccEncode(#[from] LongEccEncodeError),
        #[error(transparent)]
        LongEccDecode(#[from] LongEccDecodeError),
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

        let encoded = long::encode(&message, parity, long::OverlapFactor::Simple)?;

        // Valid data should pass validation
        assert!(validate(&encoded, parity));

        Ok(())
    }

    #[test]
    fn test_validate_long_data_invalid_with_errors() -> Result<(), TestError> {
        let message = b"This is a longer message that will use long ECC".to_buffer()?;
        let parity = 2;

        let mut encoded = long::encode(&message, parity, long::OverlapFactor::Simple)?;

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

        let encoded = long::encode(&message, parity, long::OverlapFactor::Simple)?;

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

        let encoded = long::encode(&message, parity, long::OverlapFactor::Simple)?;

        // Zero parity should still validate correctly
        assert!(validate(&encoded, parity));

        Ok(())
    }

    #[test]
    fn test_validate_long_data_header_corrupted() -> Result<(), TestError> {
        let message = b"Header corruption test".to_buffer()?;
        let parity = 2;

        let mut encoded = long::encode(&message, parity, long::OverlapFactor::Simple)?;

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
        // The message must exceed 255 bytes so that `validate` takes the long path.
        let message = b"Corrupted parity test".repeat(13).to_buffer()?;
        let parity = 2;

        let mut encoded = long::encode(&message, parity, long::OverlapFactor::Simple)?;

        // Corrupt parity bytes
        let parity_start = 32 + message.len();

        if parity_start < encoded.len() {
            encoded[parity_start] ^= 1;
        }

        // The checksum covers both the message and the parity, so corrupted
        // parity bytes invalidate the codeword.
        assert!(!validate(&encoded, parity));

        Ok(())
    }
}
