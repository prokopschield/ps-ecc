use ps_buffer::Buffer;

use crate::cow::Cow;
use crate::{Codeword, RSDecodeError, ReedSolomon};

impl ReedSolomon {
    /// Corrects a message based on detached parity bytes.
    /// # Errors
    /// - [`RSConstructorError`](crate::RSConstructorError) is returned if
    ///   `parity` holds more than [`MAX_PARITY_BYTES`](crate::MAX_PARITY_BYTES)
    ///   bytes, or an odd number of bytes.
    /// - [`std::num::TryFromIntError`] is returned if `parity` and `data`
    ///   together hold more than 255 bytes.
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSComputeErrorsError`](crate::RSComputeErrorsError) is propagated
    ///   from [`ReedSolomon::compute_errors`].
    /// - [`RSDecodeError::TooManyErrors`] is returned if the corrected bytes
    ///   fail validation.
    pub fn correct_detached<'lt>(
        parity: &[u8],
        data: &'lt [u8],
    ) -> Result<Codeword<'lt>, RSDecodeError> {
        let parity_bytes = parity.len();
        let num_parity = u8::try_from(parity_bytes >> 1)?;
        let length = u8::try_from(parity_bytes + data.len())?;
        let rs = Self::new(num_parity)?;

        let syndromes = Self::compute_syndromes_detached(parity, data)?;

        let Some(errors) = Self::compute_errors_detached(num_parity, length, &syndromes)? else {
            return Ok(data.into());
        };

        // Correct the received codeword
        let mut corrected = Buffer::with_capacity(parity.len() + data.len())?;

        corrected.extend_from_slice(parity)?;
        corrected.extend_from_slice(data)?;

        Self::apply_corrections(&mut corrected, errors.first_n_coefficients(length.into()));

        if rs.validate(&corrected).is_some() {
            return Err(RSDecodeError::TooManyErrors);
        }

        let range = parity.len()..corrected.len();
        let codeword = Cow::Owned(corrected.share());
        let codeword = Codeword { codeword, range };

        Ok(codeword)
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::{Buffer, ToBuffer};

    use crate::{RSComputeErrorsError, RSConstructorError, RSDecodeError, ReedSolomon};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_correct_detached_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DetachOk".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let corrected = ReedSolomon::correct_detached(&parity, &message)?;

        assert_eq!(&corrected[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_correct_detached_one_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"DetachErr".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut corrupted = message.clone()?;

        corrupted[1] ^= 64;

        let corrected = ReedSolomon::correct_detached(&parity, &corrupted)?;

        assert_eq!(&corrected[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_correct_detached_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"DetachMany".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut corrupted = message.clone()?;

        corrupted[0] ^= 1;
        corrupted[2] ^= 2;

        assert_eq!(
            ReedSolomon::correct_detached(&parity, &corrupted),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );

        Ok(())
    }

    #[test]
    fn test_correct_detached_rejects_oversized_parity() {
        // A 127-byte parity slice was previously accepted and silently
        // truncated to a parity count of 63; it is now rejected, matching
        // the in-place variants.
        let parity = [0u8; 127];

        assert_eq!(
            ReedSolomon::correct_detached(&parity, &[]),
            Err(RSDecodeError::RSConstructorError(
                RSConstructorError::ParityTooHigh
            ))
        );
    }

    #[test]
    fn test_correct_detached_rejects_odd_parity_length() {
        // An odd length previously dropped the last parity byte from
        // num_parity while the syndromes covered the full slice.
        let parity = [0u8; 5];

        assert_eq!(
            ReedSolomon::correct_detached(&parity, b"data"),
            Err(RSDecodeError::RSConstructorError(
                RSConstructorError::OddParityLength(5)
            ))
        );
    }

    #[test]
    fn test_correct_detached_with_errors_in_parity_only() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"ParityOnly".to_buffer()?;
        let parity_poly = rs.generate_parity(&message)?;
        let mut parity = Buffer::from_slice(parity_poly)?;

        parity[0] ^= 2;
        parity[2] ^= 4;

        let corrected = ReedSolomon::correct_detached(&parity, &message)?;

        assert_eq!(&corrected[..], &message[..]);

        Ok(())
    }
}
