use crate::finite_field::ANTILOG_TABLE;
use crate::{Polynomial, ReedSolomon};

impl ReedSolomon {
    /// Computes the syndromes of a given codeword.
    #[must_use]
    pub fn compute_syndromes(num_parity_bytes: u8, received: &[u8]) -> Polynomial {
        let num_parity_bytes = num_parity_bytes.into();

        (0..num_parity_bytes)
            .map(|i| Polynomial::eval_coefficients_at(received, ANTILOG_TABLE[i + 1].get()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::ReedSolomon;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_compute_syndromes_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Syndrome".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let syndromes = ReedSolomon::compute_syndromes(rs.parity_bytes(), &encoded);

        assert!(syndromes.is_zero());

        Ok(())
    }

    #[test]
    fn test_compute_syndromes_with_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Syndrome".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = encoded.clone()?;

        corrupted[0] ^= 1;

        let syndromes = ReedSolomon::compute_syndromes(rs.parity_bytes(), &corrupted);

        assert!(syndromes.iter().any(|&s| s != 0));

        Ok(())
    }

    #[test]
    fn test_compute_syndromes_empty_codeword() {
        let syndromes = ReedSolomon::compute_syndromes(4u8, &[]);

        assert_eq!(syndromes.degree(), 0);
        assert!(syndromes.is_zero());
    }

    /// Tests that single-error syndromes have no trailing zeros.
    ///
    /// A single error at position j produces syndrome `S_i` = e * α^(i*j), which is
    /// non-zero for all i when e ≠ 0. Thus `coefficients()` equals `first_n_coefficients()`.
    #[test]
    fn test_single_error_syndrome_length() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"Hello".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let parity_bytes: usize = rs.parity_bytes().into();

        // Single error - should always be correctable with t=3
        let mut corrupted = encoded.clone()?;

        corrupted[0] ^= 1;

        let syndromes = ReedSolomon::compute_syndromes(rs.parity_bytes(), &corrupted);

        // The full syndrome vector should have exactly parity_bytes elements
        let full = syndromes.first_n_coefficients(parity_bytes);

        assert_eq!(full.len(), parity_bytes);

        // For a single error, coefficients() has no trailing zeros to trim,
        // so it equals first_n_coefficients() in length.
        let trimmed = syndromes.coefficients();

        assert_eq!(
            trimmed.len(),
            parity_bytes,
            "coefficients() returned {} elements, expected {}",
            trimmed.len(),
            parity_bytes
        );

        Ok(())
    }
}
