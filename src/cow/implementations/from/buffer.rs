use ps_buffer::Buffer;

use crate::Cow;

impl From<Buffer> for Cow<'_> {
    fn from(value: Buffer) -> Self {
        Self::Owned(value.share())
    }
}
