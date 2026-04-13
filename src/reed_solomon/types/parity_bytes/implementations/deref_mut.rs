use std::ops::DerefMut;

use super::super::ParityBytes;

impl DerefMut for ParityBytes {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.as_mut_slice()
    }
}
