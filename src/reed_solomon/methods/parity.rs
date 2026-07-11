use crate::ReedSolomon;

impl ReedSolomon {
    /// Returns the error-correction capability of this codec: the number
    /// of byte errors a codeword can recover from.
    #[inline]
    #[must_use]
    pub const fn parity(&self) -> u8 {
        self.parity
    }
}

#[cfg(test)]
mod tests {
    use crate::{RSConstructorError, ReedSolomon};

    #[test]
    fn test_parity() -> Result<(), RSConstructorError> {
        let rs = ReedSolomon::new(8)?;

        assert_eq!(rs.parity(), 8);

        Ok(())
    }
}
