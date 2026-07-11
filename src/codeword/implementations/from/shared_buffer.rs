use ps_buffer::SharedBuffer;

use crate::Codeword;

impl From<SharedBuffer> for Codeword<'_> {
    fn from(value: SharedBuffer) -> Self {
        let range = 0..value.len();
        let codeword = value.into();

        Self { codeword, range }
    }
}
