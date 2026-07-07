/// Calculate XXH64 checksum
pub(super) fn xxh64(data: &[u8]) -> u64 {
    xxhash_rust::xxh64::xxh64(data, 8_418_112_963_040_338_442)
}
