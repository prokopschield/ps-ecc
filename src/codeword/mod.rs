mod implementations;
mod methods;

use std::ops::Range;

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
