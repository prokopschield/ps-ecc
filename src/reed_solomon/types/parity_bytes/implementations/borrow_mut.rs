use std::borrow::BorrowMut;

use super::super::ParityBytes;

impl BorrowMut<[u8]> for ParityBytes {
    fn borrow_mut(&mut self) -> &mut [u8] {
        self.as_mut()
    }
}
