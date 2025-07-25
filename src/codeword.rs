use std::ops::{Deref, Range};

use ps_buffer::Buffer;

use crate::cow::Cow;

#[derive(Debug, Hash, PartialEq, Eq)]
pub struct Codeword<'lt> {
    pub codeword: Cow<'lt>,
    pub range: Range<usize>,
}

impl<'lt> Codeword<'lt> {
    #[must_use]
    pub fn into_inner(self) -> Cow<'lt> {
        self.codeword
    }
}

impl Deref for Codeword<'_> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        &self.codeword[self.range.clone()]
    }
}

impl<'lt> From<Cow<'lt>> for Codeword<'lt> {
    fn from(value: Cow<'lt>) -> Self {
        let range = 0..value.len();
        let codeword = value;

        Self { codeword, range }
    }
}

impl<'lt> From<&'lt [u8]> for Codeword<'lt> {
    fn from(value: &'lt [u8]) -> Self {
        let range = 0..value.len();
        let codeword = value.into();

        Self { codeword, range }
    }
}

impl From<Buffer> for Codeword<'_> {
    fn from(value: Buffer) -> Self {
        let range = 0..value.len();
        let codeword = value.into();

        Self { codeword, range }
    }
}

impl AsRef<[u8]> for Codeword<'_> {
    fn as_ref(&self) -> &[u8] {
        self
    }
}
