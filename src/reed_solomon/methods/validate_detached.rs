use crate::{Polynomial, ReedSolomon};

impl ReedSolomon {
    /// Validates a segregated (parity, data) pair.
    ///
    /// Returns `None` if valid, or `Some(syndromes)` if errors are detected.
    #[must_use]
    pub fn validate_detached(parity: &[u8], data: &[u8]) -> Option<Polynomial> {
        let syndromes = Self::compute_syndromes_detached(parity, data);

        if syndromes.is_zero() {
            None
        } else {
            Some(syndromes)
        }
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::{Buffer, ToBuffer};

    use crate::ReedSolomon;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_validate_detached() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert!(ReedSolomon::validate_detached(&parity, &message).is_none());

        let mut corrupted = Buffer::with_capacity(message.len())?;

        corrupted.extend_from_slice(&message)?;
        corrupted[2] ^= 1;

        assert!(ReedSolomon::validate_detached(&parity, &corrupted).is_some());

        Ok(())
    }

    #[test]
    fn test_validate_detached_no_errors_case_2() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Detached2".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert!(ReedSolomon::validate_detached(&parity, &message).is_none());

        Ok(())
    }

    #[test]
    fn test_validate_detached_with_errors_case_2() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Detached2".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut corrupted = message.clone()?;

        corrupted[1] ^= 8;

        assert!(ReedSolomon::validate_detached(&parity, &corrupted).is_some());

        Ok(())
    }
}
