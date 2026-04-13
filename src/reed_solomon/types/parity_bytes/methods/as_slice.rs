use super::super::ParityBytes;

impl ParityBytes {
    /// Returns a slice containing the parity bytes.
    ///
    /// The returned slice has length equal to [`Self::len`].
    #[must_use]
    pub fn as_slice(&self) -> &[u8] {
        &self.data[..self.len()]
    }
}
