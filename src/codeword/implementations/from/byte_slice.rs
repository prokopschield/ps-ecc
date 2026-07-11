use crate::Codeword;

impl<'lt> From<&'lt [u8]> for Codeword<'lt> {
    fn from(value: &'lt [u8]) -> Self {
        let range = 0..value.len();
        let codeword = value.into();

        Self { codeword, range }
    }
}
