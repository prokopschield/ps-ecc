use crate::ReedSolomon;

impl ReedSolomon {
    /// XORs `corrections` onto `target`, element-wise.
    ///
    /// The two slices are zipped, so the tail of the longer one is
    /// ignored; a zero correction byte leaves its target byte unchanged.
    pub fn apply_corrections(target: &mut [u8], corrections: impl AsRef<[u8]>) {
        target
            .iter_mut()
            .zip(corrections.as_ref().iter())
            .for_each(|(target, correction)| *target ^= *correction);
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::Buffer;

    use crate::ReedSolomon;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_apply_corrections() -> Result<(), TestError> {
        let mut target = Buffer::from_slice([1, 2, 3, 4])?;
        let corrections = [0, 3, 0, 5];

        ReedSolomon::apply_corrections(&mut target, corrections);

        assert_eq!(target.as_slice(), &[1, 2 ^ 3, 3, 4 ^ 5]);

        Ok(())
    }

    #[test]
    fn test_apply_corrections_empty_target() -> Result<(), TestError> {
        let mut target = Buffer::from_slice([])?;
        let corrections = [];

        ReedSolomon::apply_corrections(&mut target, corrections);

        assert_eq!(target.as_slice(), &[]);

        Ok(())
    }
}
