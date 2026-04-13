use crate::{Polynomial, MAX_PARITY_BYTES};

use super::super::ParityBytes;

impl ParityBytes {
    /// Creates parity bytes from a polynomial remainder and parity count.
    ///
    /// Extracts the first `parity * 2` coefficients from `remainder`.
    /// The `parity` parameter is the error correction capability of the codec.
    #[must_use]
    pub fn new(remainder: &Polynomial, parity: u8) -> Self {
        let len = parity * 2;
        let coefficients = remainder.first_n_coefficients(len.into());

        let mut data = [0u8; MAX_PARITY_BYTES as usize];

        data[..coefficients.len()].copy_from_slice(coefficients);

        Self { data, len }
    }
}
