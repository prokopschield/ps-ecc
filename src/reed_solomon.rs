use ps_buffer::{Buffer, BufferError, ByteIteratorIntoBuffer, ToBuffer};

use crate::cow::Cow;
use crate::error::{
    PolynomialError, RSConstructorError, RSDecodeError, RSEncodeError, RSGenerateParityError,
};
use crate::finite_field::{div, inv, mul, ANTILOG_TABLE};
use crate::polynomial::{
    poly_div, poly_eval, poly_eval_deriv, poly_eval_detached, poly_mul, poly_rem, poly_sub,
};
use crate::{Codeword, RSComputeErrorsError, RSEuclideanError, RSValidationError};

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReedSolomon {
    parity: u8,
}

impl ReedSolomon {
    /// Creates a new Reed-Solomon codec with parameters n and k.
    /// # Errors
    /// - `ParityTooHigh` is returned if parity > 127
    pub const fn new(parity: u8) -> Result<Self, RSConstructorError> {
        use RSConstructorError::ParityTooHigh;

        if parity > 127 {
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
    /// - `PolynomialError` if generator polynomial is zero (shouldn't happen)
    pub fn generate_parity(&self, message: &[u8]) -> Result<Buffer, RSGenerateParityError> {
        let num_parity = usize::from(self.parity_bytes());
        let g = generate_generator_poly(num_parity)?;
        let dividend = vec![0u8; num_parity]
            .into_iter()
            .chain(message.iter().copied())
            .into_buffer()?;
        let mut r = poly_rem(dividend, &g)?;
        if r.len() != num_parity {
            r.resize(num_parity, 0)?;
        }
        Ok(r)
    }

    /// Encodes a message into a codeword.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - `PolynomialError` if generator polynomial is zero (shouldn't happen)
    pub fn encode(&self, message: &[u8]) -> Result<Buffer, RSEncodeError> {
        let mut buffer = Buffer::with_capacity(message.len() + usize::from(self.parity_bytes()))?;

        buffer.extend_from_slice(&self.generate_parity(message)?)?;
        buffer.extend_from_slice(message)?;

        Ok(buffer)
    }

    /// Computes the syndromes of a given codeword.
    /// # Errors
    /// - [`BufferError`] if allocation fails
    pub fn compute_syndromes(
        num_parity_bytes: impl Into<usize>,
        received: &[u8],
    ) -> Result<Buffer, BufferError> {
        (0..num_parity_bytes.into())
            .map(|i| poly_eval(received, ANTILOG_TABLE[i + 1].get()))
            .into_buffer()
    }

    /// Computes the syndromes of a given detached codeword.
    /// # Errors
    /// - [`BufferError`] if allocation fails
    pub fn compute_syndromes_detached(parity: &[u8], data: &[u8]) -> Result<Buffer, BufferError> {
        (0..parity.len())
            .map(|i| poly_eval_detached(parity, data, ANTILOG_TABLE[i + 1].get()))
            .into_buffer()
    }

    /// Validates a received codeword.
    /// # Errors
    /// `Err(syndromes)` is returned if the codeword is invalid.
    pub fn validate(&self, received: &[u8]) -> Result<Option<Buffer>, RSValidationError> {
        let syndromes = Self::compute_syndromes(self.parity_bytes(), received)?;

        if syndromes.iter().all(|&s| s == 0) {
            Ok(None)
        } else {
            Ok(Some(syndromes))
        }
    }

    /// Validates a regregated (parity, message) pair.
    /// # Errors
    /// `Err(syndromes)` is returned if the codeword is invalid.
    #[allow(clippy::cast_possible_truncation)]
    pub fn validate_detached(
        parity: &[u8],
        data: &[u8],
    ) -> Result<Option<Buffer>, RSValidationError> {
        let syndromes = Self::compute_syndromes_detached(parity, data)?;

        if syndromes.iter().all(|&s| s == 0) {
            Ok(None)
        } else {
            Ok(Some(syndromes))
        }
    }

    #[inline]
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError::TooManyErrors`] is returned if the input is unrecoverable.
    pub fn compute_errors(
        &self,
        length: usize,
        syndromes: &[u8],
    ) -> Result<Option<Buffer>, RSComputeErrorsError> {
        Self::compute_errors_detached(self.parity(), length, syndromes)
    }

    /// Computes errors in a received codeword.
    /// # Parameters
    /// - `length`: full length of codeword, including parity bytes
    /// - `syndromes`: see [`ReedSolomon::validate`]
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - `GFError` if an arithmetic operation fails
    /// - `PolynomialError` if `euclidean_for_rs` fails (division by zero)
    /// - `TooManyErrors` if the input is unrecoverable
    /// - `ZeroDerivative` shouldn't happen
    fn compute_errors_detached(
        parity: impl Into<usize>,
        length: usize,
        syndromes: &[u8],
    ) -> Result<Option<Buffer>, RSComputeErrorsError> {
        if syndromes.iter().all(|&syndrome| syndrome == 0) {
            return Ok(None);
        }

        let parity = parity.into();

        // Euclidean algorithm to find error locator and evaluator polynomials
        let (mut sigma, mut omega) = euclidean_for_rs(syndromes, parity)?;

        // Normalize sigma, omega
        let scale = inv(sigma[0])?.get();
        sigma.iter_mut().for_each(|x| *x = mul(*x, scale));
        omega.iter_mut().for_each(|x| *x = mul(*x, scale));

        // Find error positions
        let error_positions: Vec<usize> = (0..length)
            .filter(|&m| {
                let x = ANTILOG_TABLE[(255 - m) % 255].get();
                poly_eval(&sigma, x) == 0
            })
            .collect();

        if error_positions.len() > parity || error_positions.is_empty() {
            return Err(RSComputeErrorsError::TooManyErrors);
        }

        // Compute error values using Forney's formula
        let mut errors = Buffer::alloc(length)?;
        for &j in &error_positions {
            let x = ANTILOG_TABLE[(255 - j) % 255].get();
            let omega_x = poly_eval(&omega, x);
            let sigma_deriv_x = poly_eval_deriv(&sigma, x);
            if sigma_deriv_x == 0 {
                return Err(RSComputeErrorsError::ZeroErrorLocatorDerivative);
            }
            errors[j] = div(omega_x, sigma_deriv_x)?;
        }

        Ok(Some(errors))
    }

    /// Corrects a received codeword, returning the corrected codeword.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn correct<'lt>(&self, received: &'lt [u8]) -> Result<Cow<'lt>, RSDecodeError> {
        let syndromes = Self::compute_syndromes(self.parity_bytes(), received)?;

        let errors = match self.compute_errors(received.len(), &syndromes)? {
            None => return Ok(Cow::Borrowed(received)),
            Some(errors) => errors,
        };

        // Correct the received codeword
        let mut corrected = Buffer::from_slice(received)?;

        Self::apply_corrections(&mut corrected, &errors);

        match Self::validate(self, &corrected)? {
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
        let syndromes = Self::compute_syndromes(self.parity_bytes(), received)?;

        let errors = match self.compute_errors(received.len(), &syndromes)? {
            None => return Ok(()),
            Some(errors) => errors,
        };

        Self::apply_corrections(received, errors);

        match Self::validate(self, received)? {
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
        let num_parity = parity_bytes >> 1;
        let length = parity_bytes + data.len();
        let rs = Self::new(num_parity.try_into()?)?;

        let syndromes = Self::compute_syndromes_detached(parity, data)?;

        let errors = match Self::compute_errors_detached(num_parity, length, &syndromes)? {
            None => return Ok(data.into()),
            Some(errors) => errors,
        };

        // Correct the received codeword
        let mut corrected = Buffer::with_capacity(parity.len() + data.len())?;

        corrected.extend_from_slice(parity)?;
        corrected.extend_from_slice(data)?;

        Self::apply_corrections(&mut corrected, &errors);

        if let Some(_syndromes) = rs.validate(&corrected)? {
            return Err(RSDecodeError::TooManyErrors);
        };

        let range = parity.len()..corrected.len();
        let codeword = Cow::Owned(corrected);
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
        let num_parity = parity_bytes >> 1;
        let length = parity_bytes + data.len();

        let syndromes = Self::compute_syndromes_detached(parity, data)?;

        let errors = match Self::compute_errors_detached(num_parity, length, &syndromes)? {
            None => return Ok(()),
            Some(errors) => errors,
        };

        Self::apply_corrections_detached(parity, data, &errors);

        match Self::validate_detached(parity, data)? {
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

/// Generates the generator polynomial `g(x)` = `(x - α^1)(x - α^2)`...`(x - α^num_roots`).
#[allow(clippy::needless_range_loop)]
fn generate_generator_poly(num_roots: usize) -> Result<Buffer, PolynomialError> {
    let mut g = Buffer::from_slice([1])?;
    for i in 1..=num_roots {
        let root = ANTILOG_TABLE[i].get();
        g = poly_mul(&g, &[root, 1])?; // x + α^i
    }
    Ok(g)
}

/// Extended Euclidean algorithm for Reed-Solomon decoding.
fn euclidean_for_rs(s: &[u8], t: usize) -> Result<(Buffer, Buffer), RSEuclideanError> {
    let mut r0 = Buffer::alloc(2 * t + 1)?;
    r0[2 * t] = 1; // x^{2t}
    let mut r1 = s.to_buffer()?;
    let mut t0 = Buffer::from_slice([0])?;
    let mut t1 = Buffer::from_slice([1])?;
    while degree(&r1).unwrap_or(0) >= t {
        let (q, r) = poly_div(&r0, &r1)?;
        let new_t1 = poly_sub(&t0, &poly_mul(&q, &t1)?)?;
        t0 = t1;
        t1 = new_t1;
        r0 = r1;
        r1 = r;
    }
    Ok((t1, r1))
}

fn degree(poly: &[u8]) -> Option<usize> {
    poly.iter().rposition(|&x| x != 0)
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
        Polynomial(#[from] PolynomialError),
        #[error(transparent)]
        RSConstructor(#[from] RSConstructorError),
        #[error(transparent)]
        RSEncode(#[from] RSEncodeError),
        #[error(transparent)]
        RSDecode(#[from] RSDecodeError),
        #[error(transparent)]
        RSGenerateParity(#[from] RSGenerateParityError),
        #[error(transparent)]
        RSValidation(#[from] RSValidationError),
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

        assert!(rs.validate(&encoded)?.is_none());

        let mut corrupted = Buffer::with_capacity(encoded.len())?;
        corrupted.extend_from_slice(&encoded)?;
        corrupted[2] ^= 1;

        assert!(rs.validate(&corrupted)?.is_some());

        Ok(())
    }

    #[test]
    fn test_validate_detached() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Hello, World!".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert!(ReedSolomon::validate_detached(&parity, &message)?.is_none());

        let mut corrupted = Buffer::with_capacity(message.len())?;
        corrupted.extend_from_slice(&message)?;
        corrupted[2] ^= 1;

        assert!(ReedSolomon::validate_detached(&parity, &corrupted)?.is_some());

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
        assert!(ReedSolomon::new(127).is_ok());
        assert!(ReedSolomon::new(128).is_err());
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
        let syndromes = ReedSolomon::compute_syndromes(rs.parity_bytes(), &encoded)?;
        assert!(syndromes.iter().all(|&s| s == 0));
        Ok(())
    }

    #[test]
    fn test_compute_syndromes_with_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Syndrome".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[0] ^= 1;
        let syndromes = ReedSolomon::compute_syndromes(rs.parity_bytes(), &corrupted)?;
        assert!(syndromes.iter().any(|&s| s != 0));
        Ok(())
    }

    #[test]
    fn test_compute_syndromes_detached_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(3)?;
        let message = b"Detached".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &message)?;
        assert!(syndromes.iter().all(|&s| s == 0));
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
    fn test_validate_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"Valid".to_buffer()?;
        let encoded = rs.encode(&message)?;
        assert!(rs.validate(&encoded)?.is_none());
        Ok(())
    }

    #[test]
    fn test_validate_with_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(1)?;
        let message = b"Valid".to_buffer()?;
        let encoded = rs.encode(&message)?;
        let mut corrupted = encoded.clone()?;
        corrupted[0] ^= 4;
        assert!(rs.validate(&corrupted)?.is_some());
        Ok(())
    }

    #[test]
    fn test_validate_detached_no_errors_case_2() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Detached2".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        assert!(ReedSolomon::validate_detached(&parity, &message)?.is_none());
        Ok(())
    }

    #[test]
    fn test_validate_detached_with_errors_case_2() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"Detached2".to_buffer()?;
        let parity = rs.generate_parity(&message)?;
        let mut corrupted = message.clone()?;
        corrupted[1] ^= 8;
        assert!(ReedSolomon::validate_detached(&parity, &corrupted)?.is_some());
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
        assert!(rs.validate(&encoded)?.is_none());
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
    fn test_compute_syndromes_empty_codeword() -> Result<(), TestError> {
        let syndromes = ReedSolomon::compute_syndromes(4u8, &[])?;
        assert_eq!(syndromes.len(), 4);
        assert!(syndromes.iter().all(|&s| s == 0));
        Ok(())
    }

    #[test]
    fn test_compute_syndromes_detached_empty_data() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let parity = rs.generate_parity(&[])?;
        let syndromes = ReedSolomon::compute_syndromes_detached(&parity, &[])?;
        assert_eq!(syndromes.len(), 4);
        assert!(syndromes.iter().all(|&s| s == 0));
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
        assert!(rs.validate(&encoded)?.is_none());
        let mut corrupted = encoded.clone()?;
        corrupted[10] ^= 1;
        assert!(rs.validate(&corrupted)?.is_some());
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
