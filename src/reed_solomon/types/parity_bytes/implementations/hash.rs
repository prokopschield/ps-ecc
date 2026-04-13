use std::hash::{Hash, Hasher};

use super::super::ParityBytes;

impl Hash for ParityBytes {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_slice().hash(state);
    }
}
