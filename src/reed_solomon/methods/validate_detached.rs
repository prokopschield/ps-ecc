use crate::{Polynomial, RSConstructorError, ReedSolomon};

impl ReedSolomon {
    /// Validates a segregated (parity, data) pair.
    ///
    /// Returns `Ok(None)` if valid, or `Ok(Some(syndromes))` if errors are
    /// detected. An empty `parity` slice carries no parity information, so
    /// a pair within the length bound validates trivially.
    ///
    /// A pair whose combined length exceeds
    /// [`Polynomial::MAX_COEFFICIENTS`] (255) bytes can never be a valid
    /// codeword and always yields `Ok(Some)`; for such inputs the contained
    /// syndrome polynomial may be zero and is not usable for error
    /// computation.
    /// # Errors
    /// - [`RSConstructorError::ParityTooHigh`] is returned if `parity` holds
    ///   more than [`MAX_PARITY_BYTES`](crate::MAX_PARITY_BYTES) bytes.
    /// - [`RSConstructorError::OddParityLength`] is returned if `parity`
    ///   holds an odd number of bytes.
    pub fn validate_detached(
        parity: &[u8],
        data: &[u8],
    ) -> Result<Option<Polynomial>, RSConstructorError> {
        let syndromes = Self::compute_syndromes_detached(parity, data)?;

        let length_is_valid =
            parity.len() + data.len() <= usize::from(Polynomial::MAX_COEFFICIENTS);

        if length_is_valid && syndromes.is_zero() {
            Ok(None)
        } else {
            Ok(Some(syndromes))
        }
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::{Buffer, ToBuffer};

    use crate::{RSConstructorError, ReedSolomon};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_validate_detached() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert!(ReedSolomon::validate_detached(&parity, &message)?.is_none());

        let mut corrupted = Buffer::with_capacity(message.len())?;

        corrupted.extend_from_slice(&message)?;
        corrupted[2] ^= 1;

        assert!(ReedSolomon::validate_detached(&parity, &corrupted)?.is_some());

        Ok(())
    }

    #[test]
    fn test_validate_detached_no_errors_case_2() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Detached2".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert!(ReedSolomon::validate_detached(&parity, &message)?.is_none());

        Ok(())
    }

    #[test]
    fn test_validate_detached_with_errors_case_2() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Detached2".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut corrupted = message.clone()?;

        corrupted[1] ^= 8;

        assert!(ReedSolomon::validate_detached(&parity, &corrupted)?.is_some());

        Ok(())
    }

    #[test]
    fn test_validate_detached_rejects_odd_parity_length() {
        let parity = [0u8; 3];

        assert_eq!(
            ReedSolomon::validate_detached(&parity, &[]),
            Err(RSConstructorError::OddParityLength(3))
        );
    }

    #[test]
    fn test_validate_detached_empty_parity_is_trivially_valid() -> Result<(), TestError> {
        assert!(ReedSolomon::validate_detached(&[], b"anything")?.is_none());

        Ok(())
    }

    #[test]
    fn test_validate_detached_rejects_oversized_combined_length() -> Result<(), TestError> {
        // 4 parity bytes plus 252 data bytes exceed the 255-byte codeword;
        // without the upper bound the all-zero pair validated as pristine,
        // although correct_detached rejects the identical input.
        assert!(ReedSolomon::validate_detached(&[0u8; 4], &[0u8; 252])?.is_some());

        Ok(())
    }

    #[test]
    fn test_validate_detached_accepts_maximum_combined_length() -> Result<(), TestError> {
        // 4 parity bytes plus 251 data bytes fill the codeword exactly.
        assert!(ReedSolomon::validate_detached(&[0u8; 4], &[0u8; 251])?.is_none());

        Ok(())
    }

    #[test]
    fn test_validate_detached_rejects_oversized_parity() {
        // Without the length bound, an oversized parity slice was silently
        // accepted, and validation ran on silently truncated syndromes.
        for parity_len in [127usize, 256, 300] {
            let parity = vec![0u8; parity_len];

            assert_eq!(
                ReedSolomon::validate_detached(&parity, &[]),
                Err(RSConstructorError::ParityTooHigh)
            );
        }
    }
}
