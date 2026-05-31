use crate::{Polynomial, RSComputeErrorsError, ReedSolomon};

impl ReedSolomon {
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
}
