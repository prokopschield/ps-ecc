use crate::{Codeword, RSDecodeError, ReedSolomon};

impl ReedSolomon {
    /// Decodes a received codeword, correcting errors if possible.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn decode<'lt>(&self, received: &'lt [u8]) -> Result<Codeword<'lt>, RSDecodeError> {
        let corrected = self.correct(received)?;
        let codeword = Codeword {
            codeword: corrected,
            range: self.parity_bytes().into()..received.len(),
        };

        Ok(codeword)
    }
}

#[cfg(test)]
#[allow(clippy::decimal_bitwise_operands)]
mod tests {
    use ps_buffer::{Buffer, ToBuffer};

    use crate::{RSComputeErrorsError, RSDecodeError, ReedSolomon};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_encode_decode() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = Buffer::with_capacity(encoded.len())?;

        corrupted.extend_from_slice(&encoded)?;
        corrupted[2] ^= 1;

        let decoded = rs.decode(&corrupted)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Hello, World!".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = Buffer::with_capacity(encoded.len())?;

        corrupted.extend_from_slice(&encoded)?;
        corrupted[5] ^= 113;
        corrupted[6] ^= 59;
        corrupted[7] ^= 3;

        assert_eq!(
            rs.decode(&corrupted),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );

        Ok(())
    }

    #[test]
    fn test_encode_correct_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Data".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let decoded = rs.decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_correct_one_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"Example".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = encoded.clone()?;

        corrupted[1] ^= 0b1010_1010;

        let decoded = rs.decode(&corrupted)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_correct_multiple_recoverable_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Multiple".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = encoded.clone()?;

        corrupted[3] ^= 0b0011_0011;
        corrupted[7] ^= 0b1100_1100;

        let decoded = rs.decode(&corrupted)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_decode_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DecodeOk".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let decoded = rs.decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_decode_one_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"DecodeErr".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = encoded.clone()?;

        corrupted[5] ^= 32;

        let decoded = rs.decode(&corrupted)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_decode_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"DecodeMany".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = encoded.clone()?;

        corrupted[0] ^= 1;
        corrupted[2] ^= 2;

        assert_eq!(
            rs.decode(&corrupted),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );

        Ok(())
    }

    #[test]
    fn test_decode_with_errors_in_parity_only() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"ParityErr".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = encoded.clone()?;

        corrupted[1] ^= 4; // Error in parity
        corrupted[3] ^= 8; // Error in parity

        let decoded = rs.decode(&corrupted)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_correct_maximum_correctable_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"MaxErrors".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = encoded.clone()?;

        corrupted[0] ^= 1;
        corrupted[2] ^= 2;
        corrupted[4] ^= 4;
        corrupted[6] ^= 8;

        let decoded = rs.decode(&corrupted)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }
}
