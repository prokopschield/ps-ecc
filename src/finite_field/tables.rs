use std::num::NonZero;

use super::constants::{FIELD_SIZE, PRIMITIVE_POLY};

/// Log table: `log_table[x] = i` where `x = α^i`
pub const LOG_TABLE: [u8; FIELD_SIZE] = {
    let mut log = [0; FIELD_SIZE];
    let mut current = 1u16; // Start with α^0 = 1
    let mut i: u8 = 0;

    while i < 255 {
        log[current as usize] = i;
        current <<= 1; // Multiply by α (shift left in the field)
        if current & 0x100 != 0 {
            current ^= PRIMITIVE_POLY; // Reduce modulo the primitive polynomial
        }
        current &= 0xff; // Keep within 8 bits
        i += 1;
    }

    log
};

/// Declares a value as non-zero, defaults to 1
const fn nonzero(val: u8) -> NonZero<u8> {
    match NonZero::<u8>::new(val) {
        Some(value) => value,
        None => nonzero(1),
    }
}

/// Antilog table: `antilog_table[i] = α^i`
pub const ANTILOG_TABLE: [NonZero<u8>; FIELD_SIZE] = {
    let mut antilog = [nonzero(1); FIELD_SIZE];
    let mut current = 1u16; // Start with α^0 = 1
    let mut i = 0;

    while i < 255 {
        antilog[i] = nonzero((current & 0xff) as u8);
        current <<= 1;
        if current & 0x100 != 0 {
            current ^= PRIMITIVE_POLY;
        }
        i += 1;
    }

    // antilog[255] = 1; // uncomment if default value changes, see above
    // α^255 = 1 due to field order
    antilog
};

#[cfg(test)]
mod tests {
    use super::ANTILOG_TABLE;

    #[test]
    fn antilog_nonzero() {
        for n in &ANTILOG_TABLE {
            assert_ne!(n.get(), 0);
        }
    }
}
