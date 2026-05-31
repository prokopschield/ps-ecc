use crate::ReedSolomon;

impl ReedSolomon {
    #[inline]
    #[must_use]
    pub const fn parity_bytes(&self) -> u8 {
        self.parity() << 1
    }
}

#[cfg(test)]
mod tests {
    use crate::{RSConstructorError, ReedSolomon};

    #[test]
    fn test_parity_bytes() -> Result<(), RSConstructorError> {
        let rs = ReedSolomon::new(8)?;

        assert_eq!(rs.parity_bytes(), 16);

        Ok(())
    }
}
