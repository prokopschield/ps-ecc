use crate::error::GFError;

pub const FIELD_SIZE: usize = 256;
pub const PRIMITIVE_POLY: u16 = 0x11d; // x^8 + x^4 + x^3 + x^2 + 1

/// Log table: log_table[x] = i where x = α^i
pub const LOG_TABLE: [u8; FIELD_SIZE] = {
    let mut log = [0; FIELD_SIZE];
    let mut current = 1u16; // Start with α^0 = 1
    let mut i = 0;
    while i < 255 {
        log[current as usize] = i as u8;
        current <<= 1; // Multiply by α (shift left in the field)
        if current & 0x100 != 0 {
            current ^= PRIMITIVE_POLY; // Reduce modulo the primitive polynomial
        }
        current &= 0xff; // Keep within 8 bits
        i += 1;
    }
    log
};

/// Antilog table: antilog_table[i] = α^i
pub const ANTILOG_TABLE: [u8; FIELD_SIZE] = {
    let mut antilog = [0; FIELD_SIZE];
    let mut current = 1u16; // Start with α^0 = 1
    let mut i = 0;
    while i < 255 {
        antilog[i] = current as u8;
        current <<= 1;
        if current & 0x100 != 0 {
            current ^= PRIMITIVE_POLY;
        }
        current &= 0xff;
        i += 1;
    }
    antilog[255] = 1; // α^255 = 1 due to field order
    antilog
};

/// Addition in GF(256), equivalent to XOR.
pub fn add(a: u8, b: u8) -> u8 {
    a ^ b
}

/// Multiplication in GF(256) using log tables.
pub fn mul(a: u8, b: u8) -> u8 {
    if a == 0 || b == 0 {
        return 0;
    }
    let log_a = LOG_TABLE[a as usize] as u16;
    let log_b = LOG_TABLE[b as usize] as u16;
    let sum = (log_a + log_b) % 255;
    ANTILOG_TABLE[sum as usize]
}

/// Multiplicative inverse in GF(256).
pub fn inv(a: u8) -> Result<u8, GFError> {
    use GFError::*;

    if a == 0 {
        return Err(DivByZero);
    }

    let value = ANTILOG_TABLE[(255 - LOG_TABLE[a as usize] as u16) as usize];

    Ok(value)
}

/// Division in GF(256).
pub fn div(a: u8, b: u8) -> Result<u8, GFError> {
    use GFError::*;

    if b == 0 {
        return Err(DivByZero);
    }

    if a == 0 {
        return Ok(0);
    }

    Ok(mul(a, inv(b)?))
}
