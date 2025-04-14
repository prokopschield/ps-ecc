use ps_buffer::{Buffer, ToBuffer};

use crate::cow::Cow;
use crate::error::{
    PolynomialError, RSConstructorError, RSDecodeError, RSEncodeError, RSGenerateParityError,
};
use crate::finite_field::{div, inv, mul, ANTILOG_TABLE};
use crate::polynomial::{
    poly_div, poly_eval, poly_eval_deriv, poly_eval_detached, poly_mul, poly_rem, poly_sub,
};

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
    /// - `PolynomialError` if generator polynomial is zero (shouldn't happen)
    pub fn generate_parity(&self, message: &[u8]) -> Result<Vec<u8>, RSGenerateParityError> {
        let num_parity = usize::from(self.parity_bytes());
        let g = generate_generator_poly(num_parity);
        let dividend = vec![0u8; num_parity]
            .into_iter()
            .chain(message.iter().copied())
            .collect::<Vec<u8>>();
        let mut r = poly_rem(&dividend, &g)?;
        while r.len() < num_parity {
            r.push(0);
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

    /// Validates a received codeword.
    #[must_use]
    pub fn validate(&self, received: &[u8]) -> RSValidationResult {
        let num_parity_bytes = usize::from(self.parity_bytes());

        let syndromes: Vec<u8> = (0..num_parity_bytes)
            .map(|i| poly_eval(received, ANTILOG_TABLE[i + 1]))
            .collect();

        if syndromes.iter().all(|&s| s == 0) {
            RSValidationResult::Valid
        } else {
            RSValidationResult::Invalid(syndromes)
        }
    }

    /// Validates a regregated (parity, message) pair.
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn validate_detached(parity: &[u8], data: &[u8]) -> RSValidationResult {
        let syndromes: Vec<u8> = (0..parity.len())
            .map(|i| poly_eval_detached(parity, data, ANTILOG_TABLE[i + 1]))
            .collect();

        if syndromes.iter().all(|&s| s == 0) {
            RSValidationResult::Valid
        } else {
            RSValidationResult::Invalid(syndromes)
        }
    }

    #[inline]
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError::TooManyErrors`] is returned if the input is unrecoverable.
    pub fn compute_errors(&self, length: usize, syndromes: &[u8]) -> Result<Buffer, RSDecodeError> {
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
    ) -> Result<Buffer, RSDecodeError> {
        use RSDecodeError::{TooManyErrors, ZeroDerivative};

        let parity = parity.into();

        // Euclidean algorithm to find error locator and evaluator polynomials
        let (mut sigma, mut omega) = euclidean_for_rs(syndromes, parity)?;
        let scale = inv(sigma[0])?;
        sigma = sigma.iter().map(|&x| mul(x, scale)).collect();
        omega = omega.iter().map(|&x| mul(x, scale)).collect();

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
        if error_positions.len() > parity {
            return Err(TooManyErrors);
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
                return Err(ZeroDerivative);
            }
            errors[j] = div(omega_x, sigma_deriv_x)?;
        }

        Ok(errors)
    }

    /// Corrects a received codeword, returning the corrected codeword.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn correct<'lt>(&self, received: &'lt [u8]) -> Result<Cow<'lt>, RSDecodeError> {
        use RSValidationResult::{Invalid, Valid};

        // Compute syndromes
        let syndromes = match self.validate(received) {
            Valid => return Ok(Cow::Borrowed(received)),
            Invalid(syndromes) => syndromes,
        };

        let errors = self.compute_errors(received.len(), &syndromes)?;

        // Correct the received codeword
        let mut corrected = Buffer::from_slice(received)?;

        for (i, e) in errors.iter().enumerate() {
            corrected[i] ^= e;
        }

        Ok(corrected.into())
    }

    /// Decodes a received codeword, correcting errors if possible.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn decode<'lt>(&self, received: &'lt [u8]) -> Result<Cow<'lt>, RSDecodeError> {
        use RSValidationResult::{Invalid, Valid};

        let num_parity = usize::from(self.parity_bytes());

        // Compute syndromes
        let syndromes = match self.validate(received) {
            Valid => return Ok(received[num_parity..].into()),
            Invalid(syndromes) => syndromes,
        };

        let errors = self.compute_errors(received.len(), &syndromes)?;

        // Correct the received codeword
        let mut corrected = Buffer::from_slice(&received[num_parity..])?;

        for (i, e) in errors.iter().skip(num_parity).enumerate() {
            corrected[i] ^= e;
        }

        Ok(corrected.into())
    }

    /// Corrects a message based on detached parity bytes.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn correct_detached<'lt>(
        parity: &[u8],
        data: &'lt [u8],
    ) -> Result<Cow<'lt>, RSDecodeError> {
        use RSValidationResult::{Invalid, Valid};

        let parity_bytes = parity.len();
        let num_parity = parity_bytes >> 1;
        let length = parity_bytes + data.len();

        let syndromes = match Self::validate_detached(parity, data) {
            Valid => return Ok(data.into()),
            Invalid(syndromes) => syndromes,
        };

        let errors = Self::compute_errors_detached(num_parity, length, &syndromes)?;

        // Correct the received codeword
        let mut corrected = Buffer::from_slice(data)?;

        for (i, e) in errors.iter().skip(parity_bytes).enumerate() {
            corrected[i] ^= e;
        }

        Ok(corrected.into())
    }

    /// Corrects a message based on detached parity bytes.
    /// # Errors
    /// - [`RSDecodeError`] is returned if `data` is not recoverable.
    #[inline]
    pub fn correct_detached_in_place(parity: &[u8], data: &mut [u8]) -> Result<(), RSDecodeError> {
        Self::correct_both_detached_in_place(&mut parity.to_buffer()?, data)
    }

    /// Corrects a message based on detached parity bytes.
    /// Also corrects the parity bytes.
    /// # Errors
    /// - [`RSDecodeError`] is propagated from [`ReedSolomon::compute_errors`].
    pub fn correct_both_detached_in_place(
        parity: &mut [u8],
        data: &mut [u8],
    ) -> Result<(), RSDecodeError> {
        use RSValidationResult::{Invalid, Valid};

        let parity_bytes = parity.len();
        let num_parity = parity_bytes >> 1;
        let length = parity_bytes + data.len();

        let syndromes = match Self::validate_detached(parity, data) {
            Valid => return Ok(()),
            Invalid(syndromes) => syndromes,
        };

        let errors = Self::compute_errors_detached(num_parity, length, &syndromes)?;

        // Correct parity bytes
        for (i, e) in errors.iter().take(parity_bytes).enumerate() {
            parity[i] ^= e;
        }

        // Correct the detached data
        for (i, e) in errors.iter().skip(parity_bytes).enumerate() {
            data[i] ^= e;
        }

        Ok(())
    }
}

/// Generates the generator polynomial `g(x)` = `(x - α^1)(x - α^2)`...`(x - α^num_roots`).
#[allow(clippy::needless_range_loop)]
fn generate_generator_poly(num_roots: usize) -> Vec<u8> {
    let mut g = vec![1u8];
    for i in 1..=num_roots {
        let root = ANTILOG_TABLE[i];
        g = poly_mul(&g, &[root, 1]); // x + α^i
    }
    g
}

/// Extended Euclidean algorithm for Reed-Solomon decoding.
fn euclidean_for_rs(s: &[u8], t: usize) -> Result<(Vec<u8>, Vec<u8>), PolynomialError> {
    let mut r0 = vec![0u8; 2 * t + 1];
    r0[2 * t] = 1; // x^{2t}
    let mut r1 = s.to_vec();
    let mut t0 = vec![0u8];
    let mut t1 = vec![1u8];
    while degree(&r1).unwrap_or(0) >= t {
        let (q, r) = poly_div(&r0, &r1)?;
        let new_t1 = poly_sub(&t0, &poly_mul(&q, &t1));
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

pub enum RSValidationResult {
    /// Received string is correct.
    Valid,
    /// Received string is incorrect with these syndromes.
    Invalid(Vec<u8>),
}
