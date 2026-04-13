use super::super::ParityBytes;

impl Default for ParityBytes {
    fn default() -> Self {
        Self {
            data: [0u8; crate::MAX_PARITY_BYTES as usize],
            len: 0,
        }
    }
}
