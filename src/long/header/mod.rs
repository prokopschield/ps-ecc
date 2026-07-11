mod codec;
mod magic;
mod methods;
mod overlap_factor;
mod utils;

use codec::RS;
pub use magic::LONG_ECC_HEADER_MAGIC;
pub use overlap_factor::OverlapFactor;

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[allow(clippy::module_name_repetitions)]
#[repr(C, align(16))]
pub struct LongEccHeader {
    /// magic number equal to [`LONG_ECC_HEADER_MAGIC`]
    magic: u16,

    /// encoding version number (currently 1)
    version: u8,

    /// error-correction capability (each correctable error costs two parity
    /// bytes) in the low six bits, with the [`OverlapFactor`] packed into
    /// the high two bits
    parity: u8,

    /// length of the full codeword, including header and parity
    full_length: u32,

    /// length of the encoded message
    message_length: u32,

    /// XXH64 hash of the version, parity, and length fields, folded to 32 bits
    header_checksum: u32,

    /// XXH64 checksum of the message and parity bytes
    checksum: u64,

    /// parity bytes forming an RS(32, 24) codeword
    header_parity: [u8; 8],
}
