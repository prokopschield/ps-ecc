use crate::{Polynomial, ReedSolomon};

impl ReedSolomon {
    /// Validates a received codeword.
    ///
    /// Returns `None` if valid, or `Some(syndromes)` if errors are detected.
    #[must_use]
    pub fn validate(&self, received: &[u8]) -> Option<Polynomial> {
        let syndromes = Self::compute_syndromes(self.parity_bytes(), received);

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
    fn test_validate() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let encoded = rs.encode(&message)?;

        assert!(rs.validate(&encoded).is_none());

        let mut corrupted = Buffer::with_capacity(encoded.len())?;

        corrupted.extend_from_slice(&encoded)?;
        corrupted[2] ^= 1;

        assert!(rs.validate(&corrupted).is_some());

        Ok(())
    }

    #[test]
    fn test_validate_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"Valid".to_buffer()?;
        let encoded = rs.encode(&message)?;

        assert!(rs.validate(&encoded).is_none());

        Ok(())
    }

    #[test]
    fn test_validate_with_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"Valid".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = encoded.clone()?;

        corrupted[0] ^= 4;

        assert!(rs.validate(&corrupted).is_some());

        Ok(())
    }

    #[test]
    fn test_validate_large_parity() -> Result<(), TestError> {
        let rs = ReedSolomon::new(16)?;
        let message = b"LargeParity".to_buffer()?;
        let encoded = rs.encode(&message)?;

        assert!(rs.validate(&encoded).is_none());

        let mut corrupted = encoded.clone()?;

        corrupted[10] ^= 1;

        assert!(rs.validate(&corrupted).is_some());

        Ok(())
    }
}
