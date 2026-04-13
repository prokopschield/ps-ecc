use super::super::ParityBytes;

impl AsMut<[u8]> for ParityBytes {
    fn as_mut(&mut self) -> &mut [u8] {
        self.as_mut_slice()
    }
}
