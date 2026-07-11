use ps_buffer::Buffer;

use crate::cow::Cow;
use crate::{RSDecodeError, ReedSolomon};

impl ReedSolomon {
    /// Corrects a received codeword, returning the corrected codeword.
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
    pub fn correct<'lt>(&self, received: &'lt [u8]) -> Result<Cow<'lt>, RSDecodeError> {
        let parity_bytes = self.parity_bytes();

        if received.len() < usize::from(parity_bytes) {
            return Err(RSDecodeError::InsufficientLength {
                parity_bytes,
                received: received.len(),
            });
        }

        let received_len = u8::try_from(received.len())?;
        let syndromes = Self::compute_syndromes(parity_bytes, received);

        let Some(errors) = self.compute_errors(received_len, &syndromes)? else {
            return Ok(Cow::Borrowed(received));
        };

        // Correct the received codeword
        let mut corrected = Buffer::from_slice(received)?;

        Self::apply_corrections(
            &mut corrected,
            errors.first_n_coefficients(received_len.into()),
        );

        match self.validate(&corrected) {
            None => Ok(corrected.into()),
            Some(_) => Err(RSDecodeError::TooManyErrors),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{RSDecodeError, ReedSolomon};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_correct_rejects_input_shorter_than_parity() -> Result<(), TestError> {
        // An all-zero truncated slice yields zero syndromes; without the
        // length check it was returned unchanged as a pristine codeword.
        let rs = ReedSolomon::new(4)?;

        for received in [&[][..], &[0u8; 7][..]] {
            assert_eq!(
                rs.correct(received),
                Err(RSDecodeError::InsufficientLength {
                    parity_bytes: 8,
                    received: received.len(),
                })
            );
        }

        Ok(())
    }

    #[test]
    fn test_correct_accepts_parity_only_codeword() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let encoded = rs.encode(&[])?;
        let corrected = rs.correct(&encoded)?;

        assert_eq!(&corrected[..], &encoded[..]);

        Ok(())
    }
}
