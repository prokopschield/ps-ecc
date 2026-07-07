use crate::{Polynomial, RSComputeErrorsError, ReedSolomon};

impl ReedSolomon {
    /// Computes errors in a received codeword.
    /// # Parameters
    /// - `length`: full length of the codeword, including parity bytes
    /// - `syndromes`: syndrome polynomial
    /// # Errors
    /// - [`RSComputeErrorsError::GFError`] if an arithmetic operation fails
    /// - [`RSComputeErrorsError::EuclideanError`] if the Euclidean algorithm fails
    /// - [`RSComputeErrorsError::TooManyErrors`] if the input is unrecoverable
    /// - [`RSComputeErrorsError::ZeroErrorLocatorDerivative`] if the error locator derivative evaluates to zero
    #[inline]
    pub fn compute_errors(
        &self,
        length: u8,
        syndromes: &Polynomial,
    ) -> Result<Option<Polynomial>, RSComputeErrorsError> {
        Self::compute_errors_detached(self.parity(), length, syndromes)
    }
}
