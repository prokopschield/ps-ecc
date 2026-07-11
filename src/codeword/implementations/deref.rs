use std::ops::Deref;

use crate::Codeword;

impl Deref for Codeword<'_> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        &self.codeword[self.range.clone()]
    }
}
