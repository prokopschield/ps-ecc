use std::ops::Deref;

use super::super::ParityBytes;

impl Deref for ParityBytes {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.as_slice()
    }
}
