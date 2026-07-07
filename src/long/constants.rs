use super::LongEccHeader;

pub const HEADER_SIZE: usize = std::mem::size_of::<LongEccHeader>();

// Pin the wire format: the serialization methods use 32-byte arrays, so a
// layout change must fail to compile rather than panic at runtime.
const _: () = assert!(HEADER_SIZE == 32);
