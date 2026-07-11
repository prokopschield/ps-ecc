use crate::{RSDecodeError, ReedSolomon};

impl ReedSolomon {
    /// Corrects a message based on detached parity bytes.
    /// Also corrects the parity bytes.
    /// # Errors
    /// - [`RSConstructorError`](crate::RSConstructorError) is returned if
    ///   `parity` holds more than [`MAX_PARITY_BYTES`](crate::MAX_PARITY_BYTES)
    ///   bytes, or an odd number of bytes.
    /// - [`std::num::TryFromIntError`] is returned if `parity` and `data`
    ///   together hold more than 255 bytes.
    /// - [`RSComputeErrorsError`](crate::RSComputeErrorsError) is propagated
    ///   from [`ReedSolomon::compute_errors`].
    /// - [`RSDecodeError::TooManyErrors`] is returned if the corrected bytes
    ///   fail validation.
    pub fn correct_detached_in_place(
        parity: &mut [u8],
        data: &mut [u8],
    ) -> Result<(), RSDecodeError> {
        // Computing the syndromes first runs the parity-slice checks
        // (length bound, odd length) before the combined length can fail
        // the u8 conversion below.
        let syndromes = Self::compute_syndromes_detached(parity, data)?;

        let parity_bytes = parity.len();
        let num_parity = u8::try_from(parity_bytes >> 1)?;
        let length = u8::try_from(parity_bytes + data.len())?;

        let Some(errors) = Self::compute_errors_detached(num_parity, length, &syndromes)? else {
            return Ok(());
        };

        Self::apply_corrections_detached(parity, data, errors.first_n_coefficients(length.into()));

        match Self::validate_detached(parity, data)? {
            None => Ok(()),
            Some(_) => Err(RSDecodeError::TooManyErrors),
        }
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::{Buffer, ToBuffer};

    use crate::{
        RSComputeErrorsError, RSConstructorError, RSDecodeError, ReedSolomon, MAX_PARITY,
        MAX_PARITY_BYTES,
    };

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_correct_both_detached_in_place() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let mut data = Buffer::with_capacity(message.len())?;

        data.extend_from_slice(&message)?;

        let parity_poly = rs.generate_parity(&message)?;
        let mut parity = Buffer::from_slice(parity_poly)?;

        data[2] ^= 1;

        ReedSolomon::correct_detached_in_place(&mut parity, &mut data)?;

        assert_eq!(data.as_slice(), message.as_slice());
        assert_eq!(parity.as_slice(), rs.generate_parity(&message)?.as_slice());

        Ok(())
    }

    #[test]
    fn test_correct_both_detached_in_place_with_parity_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"HelloAgain!".to_buffer()?;
        let mut data = Buffer::with_capacity(message.len())?;

        data.extend_from_slice(&message)?;

        let parity_poly = rs.generate_parity(&message)?;
        let mut parity = Buffer::from_slice(parity_poly)?;

        parity[1] ^= 8;

        ReedSolomon::correct_detached_in_place(&mut parity, &mut data)?;

        assert_eq!(data.as_slice(), message.as_slice());
        assert_eq!(parity.as_slice(), rs.generate_parity(&message)?.as_slice());

        Ok(())
    }

    #[test]
    fn test_correct_both_detached_in_place_with_both_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"BothErrors".to_buffer()?;
        let mut data = Buffer::with_capacity(message.len())?;

        data.extend_from_slice(&message)?;

        let parity_poly = rs.generate_parity(&message)?;
        let mut parity = Buffer::from_slice(parity_poly)?;

        data[3] ^= 16;
        parity[0] ^= 4;

        ReedSolomon::correct_detached_in_place(&mut parity, &mut data)?;

        assert_eq!(data.as_slice(), message.as_slice());
        assert_eq!(parity.as_slice(), rs.generate_parity(&message)?.as_slice());

        Ok(())
    }

    #[test]
    fn test_correct_both_detached_in_place_max_parity() -> Result<(), TestError> {
        let rs = ReedSolomon::new(MAX_PARITY)?;
        let message = b"MaxParityInPlace".to_buffer()?;
        let mut data = message.clone()?;

        let parity_poly = rs.generate_parity(&message)?;
        let mut parity = Buffer::from_slice(parity_poly)?;

        assert_eq!(parity.len(), usize::from(MAX_PARITY_BYTES));

        data[0] ^= 0xFF;

        ReedSolomon::correct_detached_in_place(&mut parity, &mut data)?;

        assert_eq!(data.as_slice(), message.as_slice());

        Ok(())
    }

    #[test]
    fn test_correct_both_detached_in_place_rejects_oversized_parity() {
        // A parity slice with 64 nonzero bytes forms a 64-error pattern
        // relative to the all-zero codeword; without the length guard, the
        // 64th located error writes past the 63-element position array.
        for parity_len in [127usize, 128, 130] {
            let mut parity = vec![0u8; parity_len];

            parity[..64].fill(1);

            let mut data = [0u8; 10];

            let result = ReedSolomon::correct_detached_in_place(&mut parity, &mut data);

            assert_eq!(
                result,
                Err(RSDecodeError::RSConstructorError(
                    RSConstructorError::ParityTooHigh
                ))
            );
        }
    }

    #[test]
    fn test_correct_both_detached_in_place_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"TooManyBoth".to_buffer()?;
        let mut data = Buffer::with_capacity(message.len())?;

        data.extend_from_slice(&message)?;

        let parity_poly = rs.generate_parity(&message)?;
        let mut parity = Buffer::from_slice(parity_poly)?;

        data[0] ^= 1;
        data[2] ^= 2;
        parity[1] ^= 4;
        parity[3] ^= 8;

        assert_eq!(
            ReedSolomon::correct_detached_in_place(&mut parity, &mut data),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );

        Ok(())
    }
}
