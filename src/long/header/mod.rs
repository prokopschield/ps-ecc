#[derive(Clone, Copy, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[allow(clippy::module_name_repetitions)]
#[repr(C, align(16))]
pub struct LongEccHeader {
    pub full_length: u32,
    pub message_length: u32,
    pub parity: u8,
    pub segment_length: u8,
    pub segment_distance: u8,
    pub last_segment_length: u8,
    pub crc32: u32, // CRC32 of message bytes
    pub xxh64: u64, // XXH64 of message + parity
}
