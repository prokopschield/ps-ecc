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
            .collect::<Vec<u8>>();
        let mut r = poly_rem(&dividend, &g)?;
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
            .map(|i| poly_eval(received, ANTILOG_TABLE[i + 1]))
            .into_buffer()
    }

    /// Computes the syndromes of a given detached codeword.
    /// # Errors
    /// - [`BufferError`] if allocation fails
    pub fn compute_syndromes_detached(parity: &[u8], data: &[u8]) -> Result<Buffer, BufferError> {
        (0..parity.len())
            .map(|i| poly_eval_detached(parity, data, ANTILOG_TABLE[i + 1]))
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
        let scale = inv(sigma[0])?;
        sigma.iter_mut().for_each(|x| *x = mul(*x, scale));
        omega.iter_mut().for_each(|x| *x = mul(*x, scale));

        // Find error positions
        let error_positions: Vec<usize> = (0..length)
            .filter(|&m| {
                let x = if m == 0 {
                    1
                } else {
                    ANTILOG_TABLE[(255 - m) % 255]
                };
                poly_eval(&sigma, x) == 0
            })
            .collect();

        if error_positions.len() > parity || error_positions.is_empty() {
            return Err(RSComputeErrorsError::TooManyErrors);
        }

        // Compute error values using Forney's formula
        let mut errors = Buffer::alloc(length)?;
        for &j in &error_positions {
            let x = if j == 0 {
                1
            } else {
                ANTILOG_TABLE[(255 - j) % 255]
            };
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
        let root = ANTILOG_TABLE[i];
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
