use crate::{Polynomial, ReedSolomon};

impl ReedSolomon {
    /// Validates a received codeword.
    ///
    /// Returns `None` if `received` is a valid codeword, or `Some(syndromes)`
    /// if errors are detected.
    ///
    /// An input holding fewer bytes than [`ReedSolomon::parity_bytes`], or
    /// more than [`Polynomial::MAX_COEFFICIENTS`] (255) bytes, can never be
    /// a valid codeword and always yields `Some`; for such inputs the
    /// contained syndrome polynomial may be zero and is not usable for
    /// error computation.
    #[must_use]
    pub fn validate(&self, received: &[u8]) -> Option<Polynomial> {
        let parity_bytes = self.parity_bytes();
        let syndromes = Self::compute_syndromes(parity_bytes, received);

        let length_is_valid = received.len() >= usize::from(parity_bytes)
            && received.len() <= usize::from(Polynomial::MAX_COEFFICIENTS);

        if length_is_valid && syndromes.is_zero() {
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
    fn test_validate_rejects_input_shorter_than_parity() -> Result<(), TestError> {
        // An all-zero truncated slice yields zero syndromes; without the
        // length check it validated as pristine, although no codeword is
        // shorter than the parity.
        let rs = ReedSolomon::new(4)?;

        assert!(rs.validate(&[]).is_some());
        assert!(rs.validate(&[0u8; 7]).is_some());

        Ok(())
    }

    #[test]
    fn test_validate_rejects_input_longer_than_max_codeword() -> Result<(), TestError> {
        // A valid codeword zero-padded past 255 bytes yields the same zero
        // syndromes; without the upper bound it validated as pristine,
        // although decode rejects the identical bytes.
        let rs = ReedSolomon::new(4)?;

        assert!(rs.validate(&[0u8; 256]).is_some());
        assert!(rs.validate(&[0u8; 300]).is_some());

        Ok(())
    }

    #[test]
    fn test_validate_accepts_maximum_length_codeword() -> Result<(), TestError> {
        // 255 zero bytes form the codeword of the all-zero 247-byte message.
        let rs = ReedSolomon::new(4)?;

        assert!(rs.validate(&[0u8; 255]).is_none());

        Ok(())
    }

    #[test]
    fn test_validate_accepts_parity_only_codeword() -> Result<(), TestError> {
        // The encoding of the empty message is all-zero parity; a slice of
        // exactly parity_bytes zeros is that codeword and remains valid.
        let rs = ReedSolomon::new(4)?;

        assert!(rs.validate(&[0u8; 8]).is_none());

        Ok(())
    }

    #[test]
    fn test_validate_zero_parity_empty_input() -> Result<(), TestError> {
        let rs = ReedSolomon::new(0)?;

        assert!(rs.validate(&[]).is_none());

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
