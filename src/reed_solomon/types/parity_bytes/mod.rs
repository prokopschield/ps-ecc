mod implementations;
mod methods;

use crate::MAX_PARITY_BYTES;

/// Stack-allocated parity bytes returned by [`crate::ReedSolomon::generate_parity`].
#[derive(Clone, Copy, Eq)]
pub struct ParityBytes {
    /// Fixed-size buffer holding parity bytes.
    data: [u8; MAX_PARITY_BYTES as usize],

    /// Number of valid bytes in `data`.
    len: u8,
}
