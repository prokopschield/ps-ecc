/// Addition in GF(256), equivalent to XOR.
pub const fn add(a: u8, b: u8) -> u8 {
    a ^ b
}
