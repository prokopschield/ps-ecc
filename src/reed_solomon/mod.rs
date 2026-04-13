mod constants;
mod generator;

use std::ops::Rem;

pub use constants::*;
use generator::generator_poly;
use ps_buffer::{Buffer, ToBuffer};

use crate::cow::Cow;
use crate::error::{RSConstructorError, RSDecodeError, RSEncodeError, RSGenerateParityError};
use crate::finite_field::{div, inv, ANTILOG_TABLE};
use crate::{euclidean, Codeword, Polynomial, RSComputeErrorsError};

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReedSolomon {
    parity: u8,
}

impl ReedSolomon {
    /// Creates a new Reed-Solomon codec with parameters n and k.
    /// # Errors
    /// - `ParityTooHigh` is returned if parity > `MAX_PARITY`
    pub const fn new(parity: u8) -> Result<Self, RSConstructorError> {
        use RSConstructorError::ParityTooHigh;

        if parity > MAX_PARITY {
            return Err(ParityTooHigh);
        }

        let codec = Self { parity };

        Ok(codec)
    }

    #[inline]
    #[must_use]
    pub const fn parity(&self) -> u8 {
        self.parity
    }

    #[inline]
    #[must_use]
    pub const fn parity_bytes(&self) -> u8 {
        self.parity() << 1
    }

    /// Generates parity bytes.
    /// # Errors
    /// - `BufferError` if an allocation fails
    /// - [`RSGenerateParityError::Division`] if generator polynomial is zero (shouldn't happen)
    pub fn generate_parity(&self, message: &[u8]) -> Result<Buffer, RSGenerateParityError> {
        let mut p = Polynomial::default();

        p.set_coefficients(self.parity_bytes(), message)?;

        p.rem(generator_poly(self.parity()))?
            .first_n_coefficients(self.parity_bytes().into())
            .to_buffer()
            .map_err(Into::into)
    }

    /// Encodes a message into a codeword.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSEncodeError::RSGenerateParityError`] if parity generation fails.
    pub fn encode(&self, message: &[u8]) -> Result<Buffer, RSEncodeError> {
        let mut buffer = Buffer::with_capacity(message.len() + usize::from(self.parity_bytes()))?;

        buffer.extend_from_slice(&self.generate_parity(message)?)?;
        buffer.extend_from_slice(message)?;

        Ok(buffer)
    }

    /// Computes the syndromes of a given codeword.
    #[must_use]
    pub fn compute_syndromes(num_parity_bytes: u8, received: &[u8]) -> Polynomial {
        let num_parity_bytes = num_parity_bytes.into();

        (0..num_parity_bytes)
            .map(|i| Polynomial::eval_coefficients_at(received, ANTILOG_TABLE[i + 1].get()))
            .collect()
    }

    /// Computes the syndromes of a given detached codeword.
    #[must_use]
    pub fn compute_syndromes_detached(parity: &[u8], data: &[u8]) -> Polynomial {
        (0..parity.len())
            .map(|i| {
                Polynomial::eval_coefficient_slices_at(&[parity, data], ANTILOG_TABLE[i + 1].get())
            })
            .collect()
    }

    /// Validates a received codeword.
    ///
    /// Returns `None` if valid, or `Some(syndromes)` if errors are detected.
    #[must_use]
    pub fn validate(&self, received: &[u8]) -> Option<Polynomial> {
        let syndromes = Self::compute_syndromes(self.parity_bytes(), received);

        if syndromes.is_zero() {
            None
        } else {
            Some(syndromes)
        }
    }

    /// Validates a segregated (parity, data) pair.
    ///
    /// Returns `None` if valid, or `Some(syndromes)` if errors are detected.
    #[must_use]
    pub fn validate_detached(parity: &[u8], data: &[u8]) -> Option<Polynomial> {
        let syndromes = Self::compute_syndromes_detached(parity, data);

        if syndromes.is_zero() {
            None
        } else {
            Some(syndromes)
        }
    }

    /// # Errors
    /// See [`ReedSolomon::compute_errors_detached`].
    #[inline]
    pub fn compute_errors(
        &self,
        length: u8,
        syndromes: &Polynomial,
    ) -> Result<Option<Polynomial>, RSComputeErrorsError> {
        Self::compute_errors_detached(self.parity(), length, syndromes)
    }

    /// Computes errors in a received codeword.
    /// # Parameters
    /// - `length`: full length of codeword, including parity bytes
    /// - `syndromes`: syndrome polynomial
    /// # Errors
    /// - [`RSComputeErrorsError::GFError`] if an arithmetic operation fails
    /// - [`RSComputeErrorsError::EuclideanError`] if the Euclidean algorithm fails
    /// - [`RSComputeErrorsError::TooManyErrors`] if the input is unrecoverable
    /// - [`RSComputeErrorsError::ZeroErrorLocatorDerivative`] shouldn't happen
    fn compute_errors_detached(
        parity: u8,
        length: u8,
        syndromes: &Polynomial,
    ) -> Result<Option<Polynomial>, RSComputeErrorsError> {
        if syndromes.is_zero() {
            return Ok(None);
        }

        // Euclidean algorithm to find error locator and evaluator polynomials
        let (mut sigma, mut omega) = euclidean(syndromes, parity)?;

        // Normalize sigma, omega
        let scale = inv(sigma.coefficients()[0])?.get();

        sigma *= scale;
        omega *= scale;

        // Find error positions
        let mut error_positions = [0u8; MAX_PARITY as usize];
        let mut num_errors = 0usize;

        for m in 0..length {
            let x = ANTILOG_TABLE[(255 - m as usize) % 255].get();
            if Polynomial::eval_at(&sigma, x) == 0 {
                if num_errors >= usize::from(parity) {
                    return Err(RSComputeErrorsError::TooManyErrors);
                }

                error_positions[num_errors] = m;
                num_errors += 1;
            }
        }

        if num_errors == 0 {
            return Err(RSComputeErrorsError::TooManyErrors);
        }

        let error_positions = error_positions.iter().copied().take(num_errors);

        // Compute error values using Forney's formula
        let mut errors = Polynomial::default();
        for j in error_positions {
            let x = ANTILOG_TABLE[(255 - usize::from(j)) % 255].get();
            let omega_x = Polynomial::eval_at(&omega, x);
            let sigma_deriv_x = Polynomial::eval_derivative_at(&sigma, x);
            if sigma_deriv_x == 0 {
                return Err(RSComputeErrorsError::ZeroErrorLocatorDerivative);
            }
            errors.set(j, div(omega_x, sigma_deriv_x)?);
        }

        Ok(Some(errors))
    }

    /// Corrects a received codeword, returning the corrected codeword.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn correct<'lt>(&self, received: &'lt [u8]) -> Result<Cow<'lt>, RSDecodeError> {
        let received_len = u8::try_from(received.len())?;
        let syndromes = Self::compute_syndromes(self.parity_bytes(), received);

        let Some(errors) = self.compute_errors(received_len, &syndromes)? else {
            return Ok(Cow::Borrowed(received));
        };

        // Correct the received codeword
        let mut corrected = Buffer::from_slice(received)?;

        Self::apply_corrections(
            &mut corrected,
            errors.first_n_coefficients(received_len.into()),
        );

        match self.validate(&corrected) {
            None => Ok(corrected.into()),
            Some(_) => Err(RSDecodeError::TooManyErrors),
        }
    }
    /// Corrects a received codeword in-place.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    /// - [`RSDecodeError::TooManyErrors`] is returned if the data is unrecoverable.
    pub fn correct_in_place(&self, received: &mut [u8]) -> Result<(), RSDecodeError> {
        let received_len = u8::try_from(received.len())?;
        let syndromes = Self::compute_syndromes(self.parity_bytes(), received);

        let Some(errors) = self.compute_errors(received_len, &syndromes)? else {
            return Ok(());
        };

        Self::apply_corrections(received, errors.first_n_coefficients(received_len.into()));

        match self.validate(received) {
            None => Ok(()),
            Some(_) => Err(RSDecodeError::TooManyErrors),
        }
    }

    /// Decodes a received codeword, correcting errors if possible.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn decode<'lt>(&self, received: &'lt [u8]) -> Result<Codeword<'lt>, RSDecodeError> {
        let corrected = self.correct(received)?;
        let codeword = Codeword {
            codeword: corrected,
            range: self.parity_bytes().into()..received.len(),
        };

        Ok(codeword)
    }

    /// Corrects a message based on detached parity bytes.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn correct_detached<'lt>(
        parity: &[u8],
        data: &'lt [u8],
    ) -> Result<Codeword<'lt>, RSDecodeError> {
        let parity_bytes = parity.len();
        let num_parity = u8::try_from(parity_bytes >> 1)?;
        let length = u8::try_from(parity_bytes + data.len())?;
        let rs = Self::new(num_parity)?;

        let syndromes = Self::compute_syndromes_detached(parity, data);

        let Some(errors) = Self::compute_errors_detached(num_parity, length, &syndromes)? else {
            return Ok(data.into());
        };

        // Correct the received codeword
        let mut corrected = Buffer::with_capacity(parity.len() + data.len())?;

        corrected.extend_from_slice(parity)?;
        corrected.extend_from_slice(data)?;

        Self::apply_corrections(&mut corrected, errors.first_n_coefficients(length.into()));

        if rs.validate(&corrected).is_some() {
            return Err(RSDecodeError::TooManyErrors);
        }

        let range = parity.len()..corrected.len();
        let codeword = Cow::Owned(corrected.share());
        let codeword = Codeword { codeword, range };

        Ok(codeword)
    }

    /// Corrects a message based on detached parity bytes.
    /// # Errors
    /// - [`RSDecodeError`] is returned if `data` is not recoverable.
    #[inline]
    pub fn correct_detached_data_in_place(
        parity: &[u8],
        data: &mut [u8],
    ) -> Result<(), RSDecodeError> {
        Self::correct_detached_in_place(&mut parity.to_buffer()?, data)
    }

    /// Corrects a message based on detached parity bytes.
    /// Also corrects the parity bytes.
    /// # Errors
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn correct_detached_in_place(
        parity: &mut [u8],
        data: &mut [u8],
    ) -> Result<(), RSDecodeError> {
        let parity_bytes = parity.len();
        let num_parity = u8::try_from(parity_bytes >> 1)?;
        let length = u8::try_from(parity_bytes + data.len())?;

        let syndromes = Self::compute_syndromes_detached(parity, data);

        let Some(errors) = Self::compute_errors_detached(num_parity, length, &syndromes)? else {
            return Ok(());
        };

        Self::apply_corrections_detached(parity, data, errors.first_n_coefficients(length.into()));

        match Self::validate_detached(parity, data) {
            None => Ok(()),
            Some(_) => Err(RSDecodeError::TooManyErrors),
        }
    }

    pub fn apply_corrections(target: &mut [u8], corrections: impl AsRef<[u8]>) {
        target
            .iter_mut()
            .zip(corrections.as_ref().iter())
            .for_each(|(target, correction)| *target ^= *correction);
    }

    pub fn apply_corrections_detached(
        parity: &mut [u8],
        data: &mut [u8],
        corrections: impl AsRef<[u8]>,
    ) {
        let corrections = corrections.as_ref();
        Self::apply_corrections(parity, &corrections[..parity.len()]);
        Self::apply_corrections(data, &corrections[parity.len()..]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ps_buffer::ToBuffer;
    use thiserror::Error;

    #[derive(Error, Debug)]
    enum TestError {
        #[error(transparent)]
        Buffer(#[from] ps_buffer::BufferError),
        #[error(transparent)]
        RSConstructor(#[from] RSConstructorError),
        #[error(transparent)]
        RSEncode(#[from] RSEncodeError),
        #[error(transparent)]
        RSDecode(#[from] RSDecodeError),
        #[error(transparent)]
        RSGenerateParity(#[from] RSGenerateParityError),
    }

    #[test]
    fn test_encode_decode() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = Buffer::with_capacity(encoded.len())?;
        corrupted.extend_from_slice(&encoded)?;
        corrupted[2] ^= 1;

        let decoded = rs.decode(&corrupted)?;
        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Hello, World!".to_buffer()?;
        let encoded = rs.encode(&message)?;

        let mut corrupted = Buffer::with_capacity(encoded.len())?;
        corrupted.extend_from_slice(&encoded)?;
        corrupted[5] ^= 113;
        corrupted[6] ^= 59;
        corrupted[7] ^= 3;

        assert_eq!(
            rs.decode(&corrupted),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );

        Ok(())
    }

    #[test]
    fn test_validate() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let encoded = rs.encode(&message)?;

        assert!(rs.validate(&encoded).is_none());

        let mut corrupted = Buffer::with_capacity(encoded.len())?;
        corrupted.extend_from_slice(&encoded)?;
        corrupted[2] ^= 1;

        assert!(rs.validate(&corrupted).is_some());

        Ok(())
    }

    #[test]
    fn test_validate_detached() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert!(ReedSolomon::validate_detached(&parity, &message).is_none());

        let mut corrupted = Buffer::with_capacity(message.len())?;
        corrupted.extend_from_slice(&message)?;
        corrupted[2] ^= 1;

        assert!(ReedSolomon::validate_detached(&parity, &corrupted).is_some());

        Ok(())
    }

    #[test]
    fn test_correct_both_detached_in_place() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let mut data = Buffer::with_capacity(message.len())?;
        data.extend_from_slice(&message)?;
        let mut parity = rs.generate_parity(&message)?;

        data[2] ^= 1;

        ReedSolomon::correct_detached_in_place(&mut parity, &mut data)?;

        assert_eq!(data.as_slice(), message.as_slice());
        assert_eq!(parity, rs.generate_parity(&message)?);

        Ok(())
    }

    #[test]
    fn test_new() {
        assert!(ReedSolomon::new(0).is_ok());
        assert!(ReedSolomon::new(10).is_ok());
        assert!(ReedSolomon::new(MAX_PARITY).is_ok());
        assert!(ReedSolomon::new(MAX_PARITY + 1).is_err());
    }

    #[test]
    fn test_parity() -> Result<(), RSConstructorError> {
        let rs = ReedSolomon::new(8)?;
        assert_eq!(rs.parity(), 8);
        Ok(())
    }

    #[test]
    fn test_parity_bytes() -> Result<(), RSConstructorError> {
        let rs = ReedSolomon::new(8)?;
        assert_eq!(rs.parity_bytes(), 16);
        Ok(())
    }

    #[test]
    fn test_generate_parity_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Test".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        assert_eq!(parity.len(), 8); // 4 parity * 2 bytes
        Ok(())
    }

    #[test]
    fn test_generate_parity_empty_message() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        assert_eq!(parity.len(), 4);
        assert_eq!(parity.as_slice(), &[0, 0, 0, 0]);
        Ok(())
    }

    #[test]
    fn test_encode_correct_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Data".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let decoded = rs.decode(&encoded)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_encode_correct_one_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"Example".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[1] ^= 0b1010_1010;
        let decoded = rs.decode(&corrupted)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_encode_correct_multiple_recoverable_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Multiple".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[3] ^= 0b0011_0011;
        corrupted[7] ^= 0b1100_1100;
        let decoded = rs.decode(&corrupted)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

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
    fn test_validate_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"Valid".to_buffer()?;
        let encoded = rs.encode(&message)?;
        assert!(rs.validate(&encoded).is_none());
        Ok(())
    }

    #[test]
    fn test_validate_with_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"Valid".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[0] ^= 4;
        assert!(rs.validate(&corrupted).is_some());
        Ok(())
    }

    #[test]
    fn test_validate_detached_no_errors_case_2() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Detached2".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        assert!(ReedSolomon::validate_detached(&parity, &message).is_none());
        Ok(())
    }

    #[test]
    fn test_validate_detached_with_errors_case_2() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Detached2".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let mut corrupted = message.clone()?;
        corrupted[1] ^= 8;
        assert!(ReedSolomon::validate_detached(&parity, &corrupted).is_some());
        Ok(())
    }

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
    fn test_decode_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DecodeOk".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let decoded = rs.decode(&encoded)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_decode_one_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"DecodeErr".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[5] ^= 32;
        let decoded = rs.decode(&corrupted)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_decode_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"DecodeMany".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[0] ^= 1;
        corrupted[2] ^= 2;
        assert_eq!(
            rs.decode(&corrupted),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );
        Ok(())
    }

    #[test]
    fn test_correct_detached_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DetachOk".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let corrected = ReedSolomon::correct_detached(&parity, &message)?;
        assert_eq!(&corrected[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_correct_detached_one_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"DetachErr".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let mut corrupted = message.clone()?;
        corrupted[1] ^= 64;
        let corrected = ReedSolomon::correct_detached(&parity, &corrupted)?;
        assert_eq!(&corrected[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_correct_detached_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"DetachMany".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let mut corrupted = message.clone()?;
        corrupted[0] ^= 1;
        corrupted[2] ^= 2;
        assert_eq!(
            ReedSolomon::correct_detached(&parity, &corrupted),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );
        Ok(())
    }

    #[test]
    fn test_correct_detached_data_in_place_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"DataInPlaceOk".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let mut data = message.clone()?;
        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;
        assert_eq!(data.as_slice(), message.as_slice());
        Ok(())
    }

    #[test]
    fn test_correct_detached_data_in_place_one_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"DataInPlaceErr".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let mut data = message.clone()?;
        data[3] ^= 128;
        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;
        assert_eq!(data.as_slice(), message.as_slice());
        Ok(())
    }

    #[test]
    fn test_correct_detached_data_in_place_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"DataInPlaceMany".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let mut data = message.clone()?;
        data[0] ^= 1;
        data[2] ^= 2;
        assert_eq!(
            ReedSolomon::correct_detached_data_in_place(&parity, &mut data),
            Err(RSDecodeError::RSComputeErrorsError(
                RSComputeErrorsError::TooManyErrors
            ))
        );
        Ok(())
    }

    #[test]
    fn test_apply_corrections() -> Result<(), TestError> {
        let mut target = Buffer::from_slice([1, 2, 3, 4])?;
        let corrections = [0, 3, 0, 5];
        ReedSolomon::apply_corrections(&mut target, corrections);
        assert_eq!(target.as_slice(), &[1, 2 ^ 3, 3, 4 ^ 5]);
        Ok(())
    }

    #[test]
    fn test_apply_corrections_detached() -> Result<(), TestError> {
        let mut parity = Buffer::from_slice([10, 20])?;
        let mut data = Buffer::from_slice([30, 40, 50])?;
        let corrections = [1, 2, 3, 4, 5];
        ReedSolomon::apply_corrections_detached(&mut parity, &mut data, corrections);
        assert_eq!(parity.as_slice(), &[10 ^ 1, 20 ^ 2]);
        assert_eq!(data.as_slice(), &[30 ^ 3, 40 ^ 4, 50 ^ 5]);
        Ok(())
    }

    #[test]
    fn test_correct_both_detached_in_place_with_parity_error() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"HelloAgain!".to_buffer()?;
        let mut data = Buffer::with_capacity(message.len())?;
        data.extend_from_slice(&message)?;
        let mut parity = rs.generate_parity(&message)?;

        parity[1] ^= 8;

        ReedSolomon::correct_detached_in_place(&mut parity, &mut data)?;

        assert_eq!(data.as_slice(), message.as_slice());
        assert_eq!(parity, rs.generate_parity(&message)?);

        Ok(())
    }

    #[test]
    fn test_correct_both_detached_in_place_with_both_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"BothErrors".to_buffer()?;
        let mut data = Buffer::with_capacity(message.len())?;
        data.extend_from_slice(&message)?;
        let mut parity = rs.generate_parity(&message)?;

        data[3] ^= 16;
        parity[0] ^= 4;

        ReedSolomon::correct_detached_in_place(&mut parity, &mut data)?;

        assert_eq!(data.as_slice(), message.as_slice());
        assert_eq!(parity, rs.generate_parity(&message)?);

        Ok(())
    }

    #[test]
    fn test_correct_both_detached_in_place_too_many_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"TooManyBoth".to_buffer()?;
        let mut data = Buffer::with_capacity(message.len())?;
        data.extend_from_slice(&message)?;
        let mut parity = rs.generate_parity(&message)?;

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

    #[test]
    fn test_generate_parity_large_message() -> Result<(), TestError> {
        let rs = ReedSolomon::new(8)?;
        let message = vec![42; 200].to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        assert_eq!(parity.len(), 16); // 8 parity * 2 bytes
        let encoded = rs.encode(&message)?;
        assert_eq!(&encoded[..16], parity.as_slice());
        Ok(())
    }

    #[test]
    fn test_encode_empty_message() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"".to_buffer()?;
        let encoded = rs.encode(&message)?;
        assert_eq!(encoded.len(), 4); // 2 parity * 2 bytes
        assert!(rs.validate(&encoded).is_none());
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

    #[test]
    fn test_decode_with_errors_in_parity_only() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"ParityErr".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[1] ^= 4; // Error in parity
        corrupted[3] ^= 8; // Error in parity
        let decoded = rs.decode(&corrupted)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_correct_detached_data_in_place_multiple_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"MultiErr".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let mut data = message.clone()?;
        data[0] ^= 1;
        data[4] ^= 16;
        ReedSolomon::correct_detached_data_in_place(&parity, &mut data)?;
        assert_eq!(data.as_slice(), message.as_slice());
        Ok(())
    }

    #[test]
    fn test_compute_syndromes_empty_codeword() {
        let syndromes = ReedSolomon::compute_syndromes(4u8, &[]);

        assert_eq!(syndromes.degree(), 0);
        assert!(syndromes.is_zero());
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

    #[test]
    fn test_compute_syndromes_detached_empty_data() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let parity = rs.generate_parity(&[])?;
        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &[]);
        assert_eq!(syndromes.degree(), 0);
        assert!(syndromes.is_zero());
        Ok(())
    }

    #[test]
    fn test_correct_maximum_correctable_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"MaxErrors".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[0] ^= 1;
        corrupted[2] ^= 2;
        corrupted[4] ^= 4;
        corrupted[6] ^= 8;
        let decoded = rs.decode(&corrupted)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_correct_detached_with_errors_in_parity_only() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"ParityOnly".to_buffer()?;
        let mut parity = rs.generate_parity(&message)?;
        parity[0] ^= 2;
        parity[2] ^= 4;
        let corrected = ReedSolomon::correct_detached(&parity, &message)?;
        assert_eq!(&corrected[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_apply_corrections_empty_target() -> Result<(), TestError> {
        let mut target = Buffer::from_slice([])?;
        let corrections = [];
        ReedSolomon::apply_corrections(&mut target, corrections);
        assert_eq!(target.as_slice(), &[]);
        Ok(())
    }

    #[test]
    fn test_apply_corrections_detached_empty() -> Result<(), TestError> {
        let mut parity = Buffer::from_slice([])?;
        let mut data = Buffer::from_slice([])?;
        let corrections = [];
        ReedSolomon::apply_corrections_detached(&mut parity, &mut data, corrections);
        assert_eq!(parity.as_slice(), &[]);
        assert_eq!(data.as_slice(), &[]);
        Ok(())
    }

    #[test]
    fn test_validate_large_parity() -> Result<(), TestError> {
        let rs = ReedSolomon::new(16)?;
        let message = b"LargeParity".to_buffer()?;
        let encoded = rs.encode(&message)?;
        assert!(rs.validate(&encoded).is_none());
        let mut corrupted = encoded.clone()?;
        corrupted[10] ^= 1;
        assert!(rs.validate(&corrupted).is_some());
        Ok(())
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
