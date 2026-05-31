use crate::ReedSolomon;

impl ReedSolomon {
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

#[cfg(test)]
#[allow(clippy::decimal_bitwise_operands)]
mod tests {
    use ps_buffer::Buffer;

    use crate::ReedSolomon;

    type TestError = Box<dyn std::error::Error>;

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
    fn test_apply_corrections_detached_empty() -> Result<(), TestError> {
        let mut parity = Buffer::from_slice([])?;
        let mut data = Buffer::from_slice([])?;
        let corrections = [];

        ReedSolomon::apply_corrections_detached(&mut parity, &mut data, corrections);

        assert_eq!(parity.as_slice(), &[]);
        assert_eq!(data.as_slice(), &[]);

        Ok(())
    }
}
