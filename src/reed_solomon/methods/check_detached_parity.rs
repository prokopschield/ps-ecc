use crate::{RSConstructorError, ReedSolomon, MAX_PARITY_BYTES};

impl ReedSolomon {
    /// Checks that a detached parity slice is structurally valid.
    ///
    /// The checks are cheap; callers run them before any work whose cost
    /// grows with the input, such as syndrome computation.
    /// # Errors
    /// - [`RSConstructorError::ParityTooHigh`] is returned if `parity` holds
    ///   more than [`MAX_PARITY_BYTES`] bytes.
    /// - [`RSConstructorError::OddParityLength`] is returned if `parity`
    ///   holds an odd number of bytes; parity always comprises two bytes
    ///   per correctable error.
    pub(crate) fn check_detached_parity(parity: &[u8]) -> Result<(), RSConstructorError> {
        if parity.len() > usize::from(MAX_PARITY_BYTES) {
            return Err(RSConstructorError::ParityTooHigh);
        }

        if !parity.len().is_multiple_of(2) {
            return Err(RSConstructorError::OddParityLength(parity.len()));
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::{RSConstructorError, ReedSolomon, MAX_PARITY_BYTES};

    #[test]
    fn accepts_even_lengths_within_bound() {
        for parity_len in [0usize, 2, 8, MAX_PARITY_BYTES as usize] {
            let parity = vec![0u8; parity_len];

            assert_eq!(ReedSolomon::check_detached_parity(&parity), Ok(()));
        }
    }

    #[test]
    fn rejects_oversized_before_odd() {
        let parity = [0u8; 127];

        assert_eq!(
            ReedSolomon::check_detached_parity(&parity),
            Err(RSConstructorError::ParityTooHigh)
        );
    }

    #[test]
    fn rejects_odd_lengths() {
        let parity = [0u8; 5];

        assert_eq!(
            ReedSolomon::check_detached_parity(&parity),
            Err(RSConstructorError::OddParityLength(5))
        );
    }
}
