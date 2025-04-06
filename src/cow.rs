use std::ops::Deref;

use ps_buffer::{Buffer, BufferError};

pub enum Cow<'lt> {
    Borrowed(&'lt [u8]),
    Owned(Buffer),
}

impl<'lt> Cow<'lt> {
    pub fn try_into_buffer(self) -> Result<Buffer, BufferError> {
        match self {
            Cow::Borrowed(value) => Buffer::from(value),
            Cow::Owned(value) => Ok(value),
        }
    }
}

impl<'lt> Deref for Cow<'lt> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        match self {
            Self::Borrowed(value) => value,
            Self::Owned(value) => value,
        }
    }
}

impl<'lt> From<&'lt [u8]> for Cow<'lt> {
    fn from(value: &'lt [u8]) -> Self {
        Self::Borrowed(value)
    }
}

impl<'lt> From<Buffer> for Cow<'lt> {
    fn from(value: Buffer) -> Self {
        Self::Owned(value)
    }
}
