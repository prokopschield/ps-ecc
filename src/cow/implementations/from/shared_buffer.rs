use ps_buffer::SharedBuffer;

use crate::Cow;

impl From<SharedBuffer> for Cow<'_> {
    fn from(value: SharedBuffer) -> Self {
        Self::Owned(value)
    }
}
