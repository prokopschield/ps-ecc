use crate::finite_field::ANTILOG_TABLE;
use crate::{Polynomial, ReedSolomon};

impl ReedSolomon {
    /// Computes the syndromes of a given detached codeword.
    #[must_use]
    pub fn compute_syndromes_detached(parity: &[u8], data: &[u8]) -> Polynomial {
        (0..parity.len())
            .map(|i| {
                Polynomial::eval_coefficient_slices_at(&[parity, data], ANTILOG_TABLE[i + 1].get())
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::ReedSolomon;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_compute_syndromes_detached_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"Detached".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &message);

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

        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &corrupted);

        assert!(syndromes.iter().any(|&s| s != 0));

        Ok(())
    }

    #[test]
    fn test_compute_syndromes_detached_empty_data() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let parity = rs.generate_parity(&[])?;

        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &[]);

        assert_eq!(syndromes.degree(), 0);
        assert!(syndromes.is_zero());

        Ok(())
    }
}
