use super::super::ParityBytes;

impl ParityBytes {
    /// Returns `true` if the parity byte count is zero.
    ///
    /// This occurs when the Reed-Solomon codec has `parity = 0`.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }
}
