mod implementations;
mod methods;

use ps_buffer::SharedBuffer;

/// A borrowed-or-owned byte buffer.
///
/// Equality, ordering, and hashing compare the referenced bytes, so a
/// [`Cow::Borrowed`] and a [`Cow::Owned`] with the same content are equal.
#[derive(Debug)]
pub enum Cow<'lt> {
    /// Borrows the caller's bytes; no allocation took place.
    Borrowed(&'lt [u8]),
    /// Owns a shared buffer, typically holding corrected bytes.
    Owned(SharedBuffer),
}
