use crate::{RSDecodeError, ReedSolomon};

impl ReedSolomon {
    /// Corrects a received codeword in-place.
    /// # Errors
    /// - [`RSDecodeError::InsufficientLength`] is returned if `received` holds
    ///   fewer bytes than [`ReedSolomon::parity_bytes`].
    /// - [`std::num::TryFromIntError`] is returned if `received` holds more
    ///   than 255 bytes.
    /// - [`RSComputeErrorsError`](crate::RSComputeErrorsError) is propagated
    ///   from [`ReedSolomon::compute_errors`].
    /// - [`RSDecodeError::TooManyErrors`] is returned if the data is unrecoverable.
    pub fn correct_in_place(&self, received: &mut [u8]) -> Result<(), RSDecodeError> {
        let parity_bytes = self.parity_bytes();

        if received.len() < usize::from(parity_bytes) {
            return Err(RSDecodeError::InsufficientLength {
                parity_bytes,
                received: received.len(),
            });
        }

        let received_len = u8::try_from(received.len())?;
        let syndromes = Self::compute_syndromes(parity_bytes, received);

        let Some(errors) = self.compute_errors(received_len, &syndromes)? else {
            return Ok(());
        };

        Self::apply_corrections(received, errors.first_n_coefficients(received_len.into()));

        match self.validate(received) {
            None => Ok(()),
            Some(_) => Err(RSDecodeError::TooManyErrors),
        }
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::{RSComputeErrorsError, RSDecodeError, ReedSolomon};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_correct_in_place_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"InPlace".to_buffer()?;
        let mut encoded = rs.encode(&message)?;

        rs.correct_in_place(&mut encoded)?;

        assert_eq!(encoded.slice(4..), message.as_slice());

        Ok(())
    }

    #[test]
    fn test_correct_in_place_one_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"InPlace1".to_buffer()?;
        let mut encoded = rs.encode(&message)?;

        encoded[4] ^= 16;

        rs.correct_in_place(&mut encoded)?;

        assert_eq!(encoded.slice(6..), message.as_slice());

        Ok(())
    }

    #[test]
    fn test_correct_in_place_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"TooMany".to_buffer()?;
        let mut encoded = rs.encode(&message)?;

        encoded[0] ^= 1;
        encoded[1] ^= 2;

        assert_eq!(
            rs.correct_in_place(&mut encoded),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );

        Ok(())
    }

    #[test]
    fn test_correct_in_place_rejects_input_shorter_than_parity() -> Result<(), TestError> {
        // An all-zero truncated slice yields zero syndromes; without the
        // length check it was accepted unchanged as a pristine codeword.
        let rs = ReedSolomon::new(4)?;

        let mut received = [0u8; 7];

        assert_eq!(
            rs.correct_in_place(&mut received),
            Err(RSDecodeError::InsufficientLength {
                parity_bytes: 8,
                received: 7,
            })
        );

        Ok(())
    }

    #[test]
    fn test_correct_in_place_accepts_parity_only_codeword() -> Result<(), TestError> {
        // The encoding of the empty message is all-zero parity; a slice of
        // exactly parity_bytes zeros is that codeword and passes the guard.
        let rs = ReedSolomon::new(4)?;

        let mut received = [0u8; 8];

        rs.correct_in_place(&mut received)?;

        assert_eq!(received, [0u8; 8]);

        Ok(())
    }

    #[test]
    fn test_correct_in_place_multiple_errors_at_boundaries() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Boundary".to_buffer()?;
        let mut encoded = rs.encode(&message)?;
        let len = encoded.len();

        encoded[0] ^= 1; // Error at start
        encoded[len - 1] ^= 2; // Error at end

        rs.correct_in_place(&mut encoded)?;

        assert_eq!(encoded.slice(8..), message.as_slice());

        Ok(())
    }

    /// Tests that `correct_in_place` handles syndromes with trailing zeros.
    ///
    /// Although `coefficients()` trims trailing zeros, this is harmless because
    /// the euclidean algorithm converts syndromes to a `Polynomial`, which normalizes
    /// by trimming trailing zeros anyway.
    #[test]
    fn test_correct_in_place_with_trailing_zero_syndrome() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"TestMessage!".to_buffer()?;
        let encoded = rs.encode(&message)?;

        // Find a double corruption where the last syndrome is zero but others are not.
        // A single error e at position j gives syndrome S_i = e * α^(i*j), which is
        // never zero for e ≠ 0. But two errors can cancel: e1*α^(i*j1) + e2*α^(i*j2) = 0.
        'search: for pos1 in 0..encoded.len() {
            for pos2 in (pos1 + 1)..encoded.len() {
                for xor1 in 1u8..=255 {
                    for xor2 in 1u8..=255 {
                        let mut corrupted = encoded.clone()?;

                        corrupted[pos1] ^= xor1;
                        corrupted[pos2] ^= xor2;

                        let syndromes =
                            ReedSolomon::compute_syndromes(rs.parity_bytes(), &corrupted);
                        let coeffs = syndromes.first_n_coefficients(rs.parity_bytes().into());

                        // We want: last syndrome is 0, but not all are 0
                        if coeffs.last() == Some(&0) && coeffs.iter().any(|&s| s != 0) {
                            // Found a case where trailing syndrome is zero.
                            // This should still be correctable (two errors, t=2).
                            let mut to_correct = corrupted.clone()?;

                            rs.correct_in_place(&mut to_correct)?;
                            assert_eq!(to_correct.as_slice(), encoded.as_slice());
                            break 'search;
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Tests correction when syndrome polynomial has trailing zeros trimmed.
    ///
    /// Verifies that `coefficients()` being shorter than `parity_bytes` does not
    /// affect correctness, since the euclidean algorithm normalizes internally.
    #[test]
    fn test_correct_in_place_trailing_zero_syndrome_trimmed() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"TestMessage!".to_buffer()?;
        let encoded = rs.encode(&message)?;

        // Search for a corruption pattern with trailing zero syndrome
        for pos1 in 0..encoded.len() {
            for pos2 in (pos1 + 1)..encoded.len() {
                for xor1 in 1u8..=255 {
                    for xor2 in 1u8..=255 {
                        let mut corrupted = encoded.clone()?;

                        corrupted[pos1] ^= xor1;
                        corrupted[pos2] ^= xor2;

                        let syndromes =
                            ReedSolomon::compute_syndromes(rs.parity_bytes(), &corrupted);
                        let full = syndromes.first_n_coefficients(rs.parity_bytes().into());

                        if full.last() == Some(&0) && full.iter().any(|&s| s != 0) {
                            // Found pattern with trailing zero syndrome.
                            // coefficients() is shorter than parity_bytes due to trimming.
                            assert!(syndromes.coefficients().len() < rs.parity_bytes().into());

                            // Correction succeeds despite trimmed syndromes.
                            let mut to_correct = corrupted.clone()?;

                            rs.correct_in_place(&mut to_correct)?;
                            assert_eq!(to_correct.as_slice(), encoded.as_slice());

                            return Ok(());
                        }
                    }
                }
            }
        }

        panic!("No trailing zero syndrome pattern found");
    }

    #[test]
    fn test_correct_in_place_all_zeros() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = vec![0; 10].to_buffer()?;
        let mut encoded = rs.encode(&message)?;

        encoded[2] ^= 1;

        rs.correct_in_place(&mut encoded)?;

        assert_eq!(encoded.slice(4..), message.as_slice());

        Ok(())
    }
}
