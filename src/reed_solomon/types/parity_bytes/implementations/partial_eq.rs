use super::super::ParityBytes;

impl PartialEq for ParityBytes {
    fn eq(&self, other: &Self) -> bool {
        self.as_slice() == other.as_slice()
    }
}
