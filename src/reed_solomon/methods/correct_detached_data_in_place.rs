use crate::{RSConstructorError, RSDecodeError, ReedSolomon, MAX_PARITY_BYTES};

impl ReedSolomon {
    /// Corrects a message based on detached parity bytes.
    /// # Errors
    /// - [`RSConstructorError`] is returned if `parity` holds more than
    ///   [`MAX_PARITY_BYTES`] bytes, or an odd number of bytes.
    /// - [`std::num::TryFromIntError`] is returned if `parity` and `data`
    ///   together hold more than 255 bytes.
    /// - [`RSComputeErrorsError`](crate::RSComputeErrorsError) is propagated
    ///   from [`ReedSolomon::compute_errors`].
    /// - [`RSDecodeError::TooManyErrors`] is returned if the corrected bytes
    ///   fail validation.
    pub fn correct_detached_data_in_place(
        parity: &[u8],
        data: &mut [u8],
    ) -> Result<(), RSDecodeError> {
        if parity.len() > usize::from(MAX_PARITY_BYTES) {
            return Err(RSConstructorError::ParityTooHigh.into());
        }

        let parity_len = parity.len();
        let num_parity = u8::try_from(parity_len >> 1)?;
        let length = u8::try_from(parity_len + data.len())?;

        let syndromes = Self::compute_syndromes_detached(parity, data)?;

        let Some(errors) = Self::compute_errors_detached(num_parity, length, &syndromes)? else {
            return Ok(());
        };

        let corrections = errors.first_n_coefficients(length.into());

        // Apply corrections to data
        Self::apply_corrections(data, &corrections[parity_len..]);

        // Stack-allocate corrected parity for validation
        let mut corrected_parity = [0u8; MAX_PARITY_BYTES as usize];

        corrected_parity[..parity_len].copy_from_slice(parity);

        Self::apply_corrections(
            &mut corrected_parity[..parity_len],
            &corrections[..parity_len],
        );

        match Self::validate_detached(&corrected_parity[..parity_len], data)? {
            None => Ok(()),
            Some(_) => Err(RSDecodeError::TooManyErrors),
        }
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::{RSConstructorError, RSDecodeError, ReedSolomon, MAX_PARITY};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DataInPlaceOk".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn one_error_first_byte() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DataInPlaceErr".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] ^= 0x42;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn one_error_last_byte() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DataInPlaceErr".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;
        let last = data.len() - 1;

        data[last] ^= 0xFF;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn one_error_middle_byte() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DataInPlaceErr".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;
        let mid = data.len() / 2;

        data[mid] ^= 0x80;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn two_errors_max_correctable() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"TwoErrorsMax".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] ^= 0x11;
        data[5] ^= 0x22;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DataInPlaceMany".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] ^= 0x11;
        data[1] ^= 0x22;
        data[2] ^= 0x33;

        assert!(
            ReedSolomon::correct_detached_data_in_place(&parity, &mut data).is_err(),
            "Should fail with too many errors"
        );

        Ok(())
    }

    #[test]
    fn multiple_errors_higher_parity() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"HigherParityTest".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] ^= 0x01;
        data[3] ^= 0x02;
        data[7] ^= 0x04;
        data[10] ^= 0x08;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn single_byte_data() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"X".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] ^= 0xFF;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn empty_data() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn parity_error_only() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"ParityOnlyError".to_buffer()?;
        let mut parity = rs.generate_parity(&message)?;

        // Corrupt parity, not data
        parity[0] ^= 0x42;

        let mut data = message.clone()?;

        // Should still succeed because error is in parity, not data
        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn errors_in_both_parity_and_data() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"BothErrors".to_buffer()?;
        let mut parity = rs.generate_parity(&message)?;

        // Corrupt both parity and data
        parity[0] ^= 0x11;

        let mut data = message.clone()?;

        data[0] ^= 0x22;
        data[5] ^= 0x33;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn max_parity() -> Result<(), TestError> {
        let rs = ReedSolomon::new(MAX_PARITY)?;
        let message = b"MaxParityTest".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] ^= 0xFF;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn large_data() -> Result<(), TestError> {
        let rs = ReedSolomon::new(8)?;
        let message = vec![0x42u8; 200].to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] ^= 0x01;
        data[50] ^= 0x02;
        data[100] ^= 0x04;
        data[150] ^= 0x08;
        data[199] ^= 0x10;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn all_bytes_corrupted_to_zero() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"AB".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] = 0;
        data[1] = 0;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn minimum_parity() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"MinParity".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[4] ^= 0x77;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn zero_parity_no_correction() -> Result<(), TestError> {
        let rs = ReedSolomon::new(0)?;
        let message = b"ZeroParity".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert_eq!(parity.len(), 0);

        let mut data = message.clone()?;

        // With zero parity, no errors can be detected or corrected
        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn rejects_oversized_parity() {
        // Without the length guard, a 127-byte parity slice overflows the
        // 126-byte stack array used for post-correction validation, and a
        // 130-byte slice with one corrupted data byte panics while locating
        // the error past position 126.
        for parity_len in [127usize, 128, 130] {
            let parity = vec![0u8; parity_len];

            let mut data = [1u8];

            let result = ReedSolomon::correct_detached_data_in_place(&parity, &mut data);

            assert_eq!(
                result,
                Err(RSDecodeError::RSConstructorError(
                    RSConstructorError::ParityTooHigh
                ))
            );
        }
    }

    #[test]
    fn rejects_oversized_parity_with_empty_data() {
        // The combined length (130) fits in a `u8`, so only the parity
        // length guard rejects this input.
        let parity = vec![0u8; 130];

        let mut data = [];

        let result = ReedSolomon::correct_detached_data_in_place(&parity, &mut data);

        assert_eq!(
            result,
            Err(RSDecodeError::RSConstructorError(
                RSConstructorError::ParityTooHigh
            ))
        );
    }

    #[test]
    fn consecutive_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"ConsecutiveErrs".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[5] ^= 0x11;
        data[6] ^= 0x22;
        data[7] ^= 0x33;
        data[8] ^= 0x44;

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn error_values_at_boundaries() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"BoundaryVals".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        // Test with extreme error values
        let mut data = message.clone()?;

        data[0] ^= 0x01; // Minimum non-zero
        data[5] ^= 0xFF; // Maximum

        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;

        assert_eq!(&data[..], &message[..]);

        Ok(())
    }

    #[test]
    fn idempotent_correction() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Idempotent".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        let mut data = message.clone()?;

        data[0] ^= 0x42;

        // First correction
        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;
        assert_eq!(&data[..], &message[..]);

        // Second correction on already-corrected data should be no-op
        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;
        assert_eq!(&data[..], &message[..]);

        Ok(())
    }
}
