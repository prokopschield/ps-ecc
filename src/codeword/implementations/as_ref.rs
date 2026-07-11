use crate::Codeword;

impl AsRef<[u8]> for Codeword<'_> {
    fn as_ref(&self) -> &[u8] {
        self
    }
}
