use super::super::ParityBytes;

impl ParityBytes {
    /// Returns a mutable slice containing the parity bytes.
    ///
    /// The returned slice has length equal to [`Self::len`].
    #[must_use]
    pub fn as_mut_slice(&mut self) -> &mut [u8] {
        let len = self.len();
        &mut self.data[..len]
    }
}
