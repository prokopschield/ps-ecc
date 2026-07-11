/// Maximum error-correction capability of a Reed-Solomon codec.
///
/// A codec can correct at most 63 byte errors per codeword, which costs
/// 126 parity bytes: two per correctable error.
///
/// Using parity values greater than 63 is impractical in GF(256) Reed-Solomon
/// codes, as it would more than double the data size. If redundancy beyond 2x
/// is required, consider combining RS codes with replication for data storage
/// or retransmission for data transfer.
pub const MAX_PARITY: u8 = 63;

/// Maximum parity byte count, equal to `MAX_PARITY * 2`.
pub const MAX_PARITY_BYTES: u8 = MAX_PARITY * 2;
