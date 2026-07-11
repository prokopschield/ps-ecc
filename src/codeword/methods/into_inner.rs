use crate::{Codeword, Cow};

impl<'lt> Codeword<'lt> {
    /// Consumes the view and returns the full underlying buffer. Depending
    /// on the operation that produced this codeword, the buffer may carry
    /// parity and header bytes in addition to the message.
    #[must_use]
    pub fn into_inner(self) -> Cow<'lt> {
        self.codeword
    }
}
