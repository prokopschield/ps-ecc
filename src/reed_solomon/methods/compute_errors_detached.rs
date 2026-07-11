use crate::finite_field::{div, inv, ANTILOG_TABLE};
use crate::{euclidean, Polynomial, RSComputeErrorsError, ReedSolomon, MAX_PARITY};

impl ReedSolomon {
    /// Computes errors in a received codeword.
    /// # Parameters
    /// - `parity`: error-correction capability of the codec
    /// - `length`: full length of the codeword, including parity bytes
    /// - `syndromes`: syndrome polynomial
    /// # Errors
    /// - [`RSComputeErrorsError::GFError`] if an arithmetic operation fails
    /// - [`RSComputeErrorsError::EuclideanError`] if the Euclidean algorithm fails
    /// - [`RSComputeErrorsError::TooManyErrors`] if the input is unrecoverable
    /// - [`RSComputeErrorsError::ZeroErrorLocatorDerivative`] if the error locator derivative evaluates to zero
    pub(super) fn compute_errors_detached(
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
}
