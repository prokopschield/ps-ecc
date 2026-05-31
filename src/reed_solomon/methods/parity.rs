use crate::ReedSolomon;

impl ReedSolomon {
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
