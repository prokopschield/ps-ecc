use super::super::ParityBytes;

impl AsRef<[u8]> for ParityBytes {
    fn as_ref(&self) -> &[u8] {
        self.as_slice()
    }
}
