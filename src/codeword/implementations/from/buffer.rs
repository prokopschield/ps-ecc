use ps_buffer::Buffer;

use crate::Codeword;

impl From<Buffer> for Codeword<'_> {
    fn from(value: Buffer) -> Self {
        let range = 0..value.len();
        let codeword = value.into();

        Self { codeword, range }
    }
}
