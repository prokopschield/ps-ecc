use crate::{RSConstructorError, ReedSolomon, MAX_PARITY};

impl ReedSolomon {
    /// Creates a new Reed-Solomon codec with the given error-correction
    /// capability.
    /// # Errors
    /// - [`RSConstructorError::ParityTooHigh`] is returned if `parity`
    ///   exceeds [`MAX_PARITY`].
    pub const fn new(parity: u8) -> Result<Self, RSConstructorError> {
        use RSConstructorError::ParityTooHigh;

        if parity > MAX_PARITY {
            return Err(ParityTooHigh);
        }

        let codec = Self { parity };

        Ok(codec)
    }
}

#[cfg(test)]
mod tests {
    use crate::{ReedSolomon, MAX_PARITY};

    #[test]
    fn test_new() {
        assert!(ReedSolomon::new(0).is_ok());
        assert!(ReedSolomon::new(10).is_ok());
        assert!(ReedSolomon::new(MAX_PARITY).is_ok());
        assert!(ReedSolomon::new(MAX_PARITY + 1).is_err());
    }
}
