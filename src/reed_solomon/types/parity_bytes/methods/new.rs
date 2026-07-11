use crate::{Polynomial, MAX_PARITY, MAX_PARITY_BYTES};

use super::super::ParityBytes;

impl ParityBytes {
    /// Creates parity bytes from a polynomial remainder and parity count.
    ///
    /// Extracts the first `parity * 2` coefficients from `remainder`.
    /// The `parity` parameter is the error correction capability of the codec
    /// and must not exceed [`MAX_PARITY`]; the sole caller,
    /// [`ReedSolomon::generate_parity`](crate::ReedSolomon::generate_parity),
    /// passes a parity count already validated by
    /// [`ReedSolomon::new`](crate::ReedSolomon::new).
    /// # Panics
    /// Panics in debug builds if `parity` exceeds [`MAX_PARITY`].
    #[must_use]
    pub(crate) fn new(remainder: &Polynomial, parity: u8) -> Self {
        debug_assert!(
            parity <= MAX_PARITY,
            "parity {parity} exceeds MAX_PARITY {MAX_PARITY}"
        );

        let len = parity * 2;
        let coefficients = remainder.first_n_coefficients(len.into());

        let mut data = [0u8; MAX_PARITY_BYTES as usize];

        data[..coefficients.len()].copy_from_slice(coefficients);

        Self { data, len }
    }
}

#[cfg(test)]
mod tests {
    use crate::{ParityBytes, Polynomial, MAX_PARITY, MAX_PARITY_BYTES};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn copies_two_coefficients_per_parity_symbol() -> Result<(), TestError> {
        let remainder = Polynomial::from_slice(&[1, 2, 3, 4, 5])?;
        let parity_bytes = ParityBytes::new(&remainder, 2);

        assert_eq!(parity_bytes.len(), 4);
        assert_eq!(parity_bytes.as_slice(), &[1, 2, 3, 4]);

        Ok(())
    }

    #[test]
    fn zero_parity_is_empty() {
        let remainder = Polynomial::default();
        let parity_bytes = ParityBytes::new(&remainder, 0);

        assert!(parity_bytes.is_empty());
    }

    #[test]
    fn max_parity_fills_entire_buffer() -> Result<(), TestError> {
        let remainder = Polynomial::from_slice(&[7u8; MAX_PARITY_BYTES as usize])?;
        let parity_bytes = ParityBytes::new(&remainder, MAX_PARITY);

        assert_eq!(parity_bytes.len(), usize::from(MAX_PARITY_BYTES));
        assert!(parity_bytes.as_slice().iter().all(|&b| b == 7));

        Ok(())
    }

    #[test]
    #[should_panic(expected = "exceeds MAX_PARITY")]
    #[cfg(debug_assertions)]
    fn rejects_parity_above_max_in_debug_builds() {
        let remainder = Polynomial::default();
        let _ = ParityBytes::new(&remainder, MAX_PARITY + 1);
    }
}
