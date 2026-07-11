use std::ops::{Deref, Range};

use ps_buffer::{Buffer, SharedBuffer};

use crate::cow::Cow;

/// A decoded codeword exposing the message bytes.
///
/// Dereferences to the message portion of the codeword, selected by
/// `range`. The fields are crate-internal so that a range outside the
/// codeword, which would make dereferencing panic, cannot be constructed
/// by callers.
#[derive(Debug, Hash, PartialEq, Eq)]
pub struct Codeword<'lt> {
    pub(crate) codeword: Cow<'lt>,
    pub(crate) range: Range<usize>,
}

impl<'lt> Codeword<'lt> {
    /// Consumes the view and returns the full underlying codeword,
    /// including the parity bytes and, for long codewords, the header.
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

impl From<SharedBuffer> for Codeword<'_> {
    fn from(value: SharedBuffer) -> Self {
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
