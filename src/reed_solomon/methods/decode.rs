use crate::{Codeword, RSDecodeError, ReedSolomon};

impl ReedSolomon {
    /// Decodes a received codeword, correcting errors if possible.
    /// # Errors
    /// - [`RSDecodeError::InsufficientLength`] is returned if `received` holds
    ///   fewer bytes than [`ReedSolomon::parity_bytes`].
    /// - [`std::num::TryFromIntError`] is returned if `received` holds more
    ///   than 255 bytes.
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSComputeErrorsError`](crate::RSComputeErrorsError) is propagated
    ///   from [`ReedSolomon::compute_errors`].
    /// - [`RSDecodeError::TooManyErrors`] is returned if the corrected bytes
    ///   fail validation.
    pub fn decode<'lt>(&self, received: &'lt [u8]) -> Result<Codeword<'lt>, RSDecodeError> {
        let parity_bytes = self.parity_bytes();

        if received.len() < usize::from(parity_bytes) {
            return Err(RSDecodeError::InsufficientLength {
                parity_bytes,
                received: received.len(),
            });
        }

        let corrected = self.correct(received)?;
        let codeword = Codeword {
            codeword: corrected,
            range: usize::from(parity_bytes)..received.len(),
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
    fn test_decode_empty_input_rejected() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;

        // Without the length guard, decode returns Ok with the inverted
        // range 8..0, and the first deref panics.
        assert_eq!(
            rs.decode(&[]),
            Err(RSDecodeError::InsufficientLength {
                parity_bytes: 8,
                received: 0,
            })
        );

        Ok(())
    }

    #[test]
    fn test_decode_input_shorter_than_parity_rejected() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;

        assert_eq!(
            rs.decode(&[0u8; 7]),
            Err(RSDecodeError::InsufficientLength {
                parity_bytes: 8,
                received: 7,
            })
        );

        Ok(())
    }

    #[test]
    fn test_decode_parity_only_codeword() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let encoded = rs.encode(&[])?;

        assert_eq!(encoded.len(), usize::from(rs.parity_bytes()));

        let decoded = rs.decode(&encoded)?;

        assert!(decoded.is_empty());

        Ok(())
    }

    #[test]
    fn test_decode_zero_parity_empty_input() -> Result<(), TestError> {
        let rs = ReedSolomon::new(0)?;
        let decoded = rs.decode(&[])?;

        assert!(decoded.is_empty());

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
