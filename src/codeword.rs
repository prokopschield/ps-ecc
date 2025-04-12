use std::ops::{Deref, Range};

use ps_buffer::Buffer;

use crate::cow::Cow;

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

impl<'lt> Deref for Codeword<'lt> {
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

impl<'lt> From<Buffer> for Codeword<'lt> {
    fn from(value: Buffer) -> Self {
        let range = 0..value.len();
        let codeword = value.into();

        Self { codeword, range }
    }
}
