use crate::{Polynomial, RSConstructorError, ReedSolomon};

impl ReedSolomon {
    /// Validates a segregated (parity, data) pair.
    ///
    /// Returns `Ok(None)` if valid, or `Ok(Some(syndromes))` if errors are
    /// detected.
    /// # Errors
    /// - [`RSConstructorError::ParityTooHigh`] is returned if `parity` holds
    ///   more than [`MAX_PARITY_BYTES`](crate::MAX_PARITY_BYTES) bytes.
    pub fn validate_detached(
        parity: &[u8],
        data: &[u8],
    ) -> Result<Option<Polynomial>, RSConstructorError> {
        let syndromes = Self::compute_syndromes_detached(parity, data)?;

        if syndromes.is_zero() {
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
