use crate::finite_field::ANTILOG_TABLE;
use crate::{Polynomial, RSConstructorError, ReedSolomon, MAX_PARITY_BYTES};

impl ReedSolomon {
    /// Computes the syndromes of a given detached codeword.
    ///
    /// An empty `parity` slice yields zero syndromes, so a pair carrying no
    /// parity always validates as pristine.
    /// # Errors
    /// - [`RSConstructorError::ParityTooHigh`] is returned if `parity` holds
    ///   more than [`MAX_PARITY_BYTES`] bytes.
    /// - [`RSConstructorError::OddParityLength`] is returned if `parity`
    ///   holds an odd number of bytes; parity is always generated in
    ///   two-byte symbols.
    pub fn compute_syndromes_detached(
        parity: &[u8],
        data: &[u8],
    ) -> Result<Polynomial, RSConstructorError> {
        if parity.len() > usize::from(MAX_PARITY_BYTES) {
            return Err(RSConstructorError::ParityTooHigh);
        }

        if !parity.len().is_multiple_of(2) {
            return Err(RSConstructorError::OddParityLength(parity.len()));
        }

        let syndromes = (0..parity.len())
            .map(|i| {
                Polynomial::eval_coefficient_slices_at(&[parity, data], ANTILOG_TABLE[i + 1].get())
            })
            .collect();

        Ok(syndromes)
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::{RSConstructorError, ReedSolomon, MAX_PARITY_BYTES};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_compute_syndromes_detached_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"Detached".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &message)?;

        assert!(syndromes.is_zero());

        Ok(())
    }

    #[test]
    fn test_compute_syndromes_detached_with_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"Detached".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut corrupted = message.clone()?;

        corrupted[2] ^= 2;

        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &corrupted)?;

        assert!(syndromes.iter().any(|&s| s != 0));

        Ok(())
    }

    #[test]
    fn test_compute_syndromes_detached_empty_data() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let parity = rs.generate_parity(&[])?;

        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &[])?;

        assert_eq!(syndromes.degree(), 0);
        assert!(syndromes.is_zero());

        Ok(())
    }

    #[test]
    fn test_compute_syndromes_detached_max_parity_bytes() -> Result<(), TestError> {
        let parity = [0u8; MAX_PARITY_BYTES as usize];

        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &[])?;

        assert!(syndromes.is_zero());

        Ok(())
    }

    #[test]
    fn test_compute_syndromes_detached_rejects_odd_parity_length() {
        // num_parity = parity.len() >> 1 silently dropped the odd byte
        // while the syndromes covered the full slice; odd lengths are now
        // rejected because generate_parity only produces two-byte symbols.
        for parity_len in [1usize, 3, 125] {
            let parity = vec![0u8; parity_len];

            assert_eq!(
                ReedSolomon::compute_syndromes_detached(&parity, &[]),
                Err(RSConstructorError::OddParityLength(parity_len))
            );
        }
    }

    #[test]
    fn test_compute_syndromes_detached_oversized_wins_over_odd() {
        // 127 is both oversized and odd; the length bound is checked first.
        let parity = [0u8; 127];

        assert_eq!(
            ReedSolomon::compute_syndromes_detached(&parity, &[]),
            Err(RSConstructorError::ParityTooHigh)
        );
    }

    #[test]
    fn test_compute_syndromes_detached_rejects_oversized_parity() {
        // Without the length bound, an oversized parity slice was silently
        // accepted; collection into Polynomial consumes at most 255
        // coefficients, so the syndromes were silently truncated.
        for parity_len in [127usize, 256, 300] {
            let parity = vec![0u8; parity_len];

            assert_eq!(
                ReedSolomon::compute_syndromes_detached(&parity, &[]),
                Err(RSConstructorError::ParityTooHigh)
            );
        }
    }
}
