use crate::error::{PolynomialError, RSConstructorError, RSDecodeError, RSEncodeError};
use crate::finite_field::{add, div, inv, mul, ANTILOG_TABLE};
use crate::polynomial::{poly_div, poly_eval, poly_eval_deriv, poly_mul, poly_rem};

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
    pub fn generate_parity(&self, message: &[u8]) -> Result<Vec<u8>, PolynomialError> {
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
    /// - `PolynomialError` if generator polynomial is zero (shouldn't happen)
    pub fn encode(&self, message: &[u8]) -> Result<Vec<u8>, RSEncodeError> {
        let mut r = self.generate_parity(message)?;
        r.extend_from_slice(message);
        Ok(r)
    }

    /// Decodes a received codeword, correcting errors if possible.
    /// # Errors
    /// - `GFError` if an arithmetic operation fails
    /// - `PolynomialError` if `euclidean_for_rs` fails (division by zero)
    /// - `TooManyErrors` if the input is unrecoverable
    /// - `ZeroDerivative` shouldn't happen
    pub fn decode(&self, received: &[u8]) -> Result<Vec<u8>, RSDecodeError> {
        use RSDecodeError::{TooManyErrors, ZeroDerivative};

        let num_parity = usize::from(self.parity_bytes());
        let t = usize::from(self.parity());
        let n = received.len();

        // Compute syndromes
        let mut syndromes = vec![0u8; num_parity];
        for i in 0..num_parity {
            syndromes[i] = poly_eval(received, ANTILOG_TABLE[i + 1]);
        }
        if syndromes.iter().all(|&s| s == 0) {
            return Ok(received[num_parity..].to_vec());
        }

        // Euclidean algorithm to find error locator and evaluator polynomials
        let (mut sigma, mut omega) = euclidean_for_rs(&syndromes, t)?;
        let scale = inv(sigma[0])?;
        sigma = sigma.iter().map(|&x| mul(x, scale)).collect();
        omega = omega.iter().map(|&x| mul(x, scale)).collect();

        // Find error positions
        let error_positions: Vec<usize> = (0..n)
            .filter(|&m| {
                let x = if m == 0 {
                    1
                } else {
                    ANTILOG_TABLE[(255 - m) % 255]
                };
                poly_eval(&sigma, x) == 0
            })
            .collect();
        if error_positions.len() > t {
            return Err(TooManyErrors);
        }

        // Compute error values using Forney's formula
        let mut errors = vec![0u8; n];
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

        // Correct the received codeword
        let corrected = received
            .iter()
            .zip(errors.iter())
            .map(|(&r, &e)| add(r, e))
            .collect::<Vec<u8>>();
        Ok(corrected[num_parity..].to_vec())
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

/// Subtracts two polynomials (same as addition in GF(2)).
#[allow(clippy::needless_range_loop)]
fn poly_sub(p1: &[u8], p2: &[u8]) -> Vec<u8> {
    let len = p1.len().max(p2.len());
    let mut result = vec![0u8; len];
    for i in 0..len {
        let a = p1.get(i).copied().unwrap_or(0);
        let b = p2.get(i).copied().unwrap_or(0);
        result[i] = add(a, b);
    }
    trim_leading_zeros(&mut result);
    result
}

fn degree(poly: &[u8]) -> Option<usize> {
    poly.iter().rposition(|&x| x != 0)
}

fn trim_leading_zeros(poly: &mut Vec<u8>) {
    while poly.len() > 1 && poly.last() == Some(&0) {
        poly.pop();
    }
}
