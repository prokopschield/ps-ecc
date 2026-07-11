use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::ops::Deref;

use ps_buffer::{Buffer, BufferError, SharedBuffer};

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

impl PartialEq for Cow<'_> {
    fn eq(&self, other: &Self) -> bool {
        **self == **other
    }
}

impl Eq for Cow<'_> {}

impl PartialOrd for Cow<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cow<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        (**self).cmp(&**other)
    }
}

impl Hash for Cow<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (**self).hash(state);
    }
}

impl Cow<'_> {
    /// Converts this [`Cow`] into a [`SharedBuffer`].
    /// - In the case of [`Cow::Borrowed`], a new buffer is allocated.
    /// - In the case of [`Cow::Owned`], the existing buffer is returned.
    /// # Errors
    /// [`BufferError`] is returned if an allocation error occurs.
    pub fn try_into_buffer(self) -> Result<SharedBuffer, BufferError> {
        match self {
            Cow::Borrowed(value) => Ok(Buffer::from_slice(value)?.share()),
            Cow::Owned(value) => Ok(value),
        }
    }
}

impl Deref for Cow<'_> {
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

impl From<Buffer> for Cow<'_> {
    fn from(value: Buffer) -> Self {
        Self::Owned(value.share())
    }
}

impl From<SharedBuffer> for Cow<'_> {
    fn from(value: SharedBuffer) -> Self {
        Self::Owned(value)
    }
}

#[cfg(test)]
mod tests {
    use std::hash::{DefaultHasher, Hash, Hasher};

    use ps_buffer::Buffer;

    use super::Cow;

    type TestError = Box<dyn std::error::Error>;

    fn hash(value: &Cow) -> u64 {
        let mut hasher = DefaultHasher::new();

        value.hash(&mut hasher);

        hasher.finish()
    }

    #[test]
    fn test_eq_across_variants_with_same_content() -> Result<(), TestError> {
        let borrowed = Cow::Borrowed(b"abc");
        let owned = Cow::Owned(Buffer::from_slice(b"abc")?.share());

        assert_eq!(borrowed, owned);

        Ok(())
    }

    #[test]
    fn test_ne_for_different_content() -> Result<(), TestError> {
        let borrowed = Cow::Borrowed(b"abc");
        let owned = Cow::Owned(Buffer::from_slice(b"abd")?.share());

        assert_ne!(borrowed, owned);

        Ok(())
    }

    #[test]
    fn test_hash_agrees_across_variants() -> Result<(), TestError> {
        let borrowed = Cow::Borrowed(b"abc");
        let owned = Cow::Owned(Buffer::from_slice(b"abc")?.share());

        assert_eq!(hash(&borrowed), hash(&owned));

        Ok(())
    }

    #[test]
    fn test_ord_compares_content_across_variants() -> Result<(), TestError> {
        let borrowed = Cow::Borrowed(b"b");
        let owned = Cow::Owned(Buffer::from_slice(b"a")?.share());

        assert!(borrowed > owned);
        assert!(owned < borrowed);

        Ok(())
    }
}
