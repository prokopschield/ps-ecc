use std::borrow::Borrow;

use super::super::ParityBytes;

impl Borrow<[u8]> for ParityBytes {
    fn borrow(&self) -> &[u8] {
        self.as_slice()
    }
}
