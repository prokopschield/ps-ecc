use ps_buffer::Buffer;

use crate::cow::Cow;
use crate::{RSDecodeError, ReedSolomon};

impl ReedSolomon {
    /// Corrects a received codeword, returning the corrected codeword.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn correct<'lt>(&self, received: &'lt [u8]) -> Result<Cow<'lt>, RSDecodeError> {
        let received_len = u8::try_from(received.len())?;
        let syndromes = Self::compute_syndromes(self.parity_bytes(), received);

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
