use ps_buffer::{Buffer, BufferError, SharedBuffer};

use crate::Cow;

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
