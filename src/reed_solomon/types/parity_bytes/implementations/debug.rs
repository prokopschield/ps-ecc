use std::fmt::{Debug, Formatter, Result};

use super::super::ParityBytes;

impl Debug for ParityBytes {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        self.as_slice().fmt(f)
    }
}
