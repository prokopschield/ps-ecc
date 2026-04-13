use super::super::ParityBytes;

impl ParityBytes {
    /// Returns the number of parity bytes.
    ///
    /// Equal to `parity * 2`, where `parity` is the error correction capability.
    /// The value ranges from 0 to [`crate::MAX_PARITY_BYTES`].
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len as usize
    }
}
