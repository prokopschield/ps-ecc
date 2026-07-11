use crate::ReedSolomon;

impl ReedSolomon {
    /// Applies XOR corrections to a detached `(parity, data)` pair.
    ///
    /// The `corrections` slice is laid out parity-first, matching the
    /// `parity || data` codeword layout: its first `parity.len()` bytes are
    /// `XORed` onto `parity`, and the remaining bytes onto `data`.
    ///
    /// Like [`ReedSolomon::apply_corrections`], this method uses truncating
    /// zip semantics. If `corrections` holds fewer bytes than `parity` and
    /// `data` combined, the uncovered tail is left unmodified; excess
    /// correction bytes are ignored.
    pub fn apply_corrections_detached(
        parity: &mut [u8],
        data: &mut [u8],
        corrections: impl AsRef<[u8]>,
    ) {
        let corrections = corrections.as_ref();
        let split = corrections.len().min(parity.len());
        let (parity_corrections, data_corrections) = corrections.split_at(split);

        Self::apply_corrections(parity, parity_corrections);
        Self::apply_corrections(data, data_corrections);
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

    #[test]
    fn test_corrections_shorter_than_parity() -> Result<(), TestError> {
        // Previously panicked on `corrections[..parity.len()]`; now the
        // uncovered parity tail and all of `data` are left unmodified.
        let mut parity = Buffer::from_slice([10, 20, 30])?;
        let mut data = Buffer::from_slice([40, 50])?;
        let corrections = [1];

        ReedSolomon::apply_corrections_detached(&mut parity, &mut data, corrections);

        assert_eq!(parity.as_slice(), &[10 ^ 1, 20, 30]);
        assert_eq!(data.as_slice(), &[40, 50]);

        Ok(())
    }

    #[test]
    fn test_corrections_cover_parity_only() -> Result<(), TestError> {
        let mut parity = Buffer::from_slice([10, 20])?;
        let mut data = Buffer::from_slice([30, 40])?;
        let corrections = [1, 2];

        ReedSolomon::apply_corrections_detached(&mut parity, &mut data, corrections);

        assert_eq!(parity.as_slice(), &[10 ^ 1, 20 ^ 2]);
        assert_eq!(data.as_slice(), &[30, 40]);

        Ok(())
    }

    #[test]
    fn test_corrections_cover_parity_and_partial_data() -> Result<(), TestError> {
        let mut parity = Buffer::from_slice([10, 20])?;
        let mut data = Buffer::from_slice([30, 40, 50])?;
        let corrections = [1, 2, 3];

        ReedSolomon::apply_corrections_detached(&mut parity, &mut data, corrections);

        assert_eq!(parity.as_slice(), &[10 ^ 1, 20 ^ 2]);
        assert_eq!(data.as_slice(), &[30 ^ 3, 40, 50]);

        Ok(())
    }

    #[test]
    fn test_corrections_longer_than_codeword() -> Result<(), TestError> {
        let mut parity = Buffer::from_slice([10])?;
        let mut data = Buffer::from_slice([20])?;
        let corrections = [1, 2, 3, 4];

        ReedSolomon::apply_corrections_detached(&mut parity, &mut data, corrections);

        assert_eq!(parity.as_slice(), &[10 ^ 1]);
        assert_eq!(data.as_slice(), &[20 ^ 2]);

        Ok(())
    }

    #[test]
    fn test_empty_parity_with_corrections() -> Result<(), TestError> {
        let mut parity = Buffer::from_slice([])?;
        let mut data = Buffer::from_slice([30, 40])?;
        let corrections = [1, 2];

        ReedSolomon::apply_corrections_detached(&mut parity, &mut data, corrections);

        assert_eq!(parity.as_slice(), &[]);
        assert_eq!(data.as_slice(), &[30 ^ 1, 40 ^ 2]);

        Ok(())
    }
}
