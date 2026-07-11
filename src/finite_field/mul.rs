use super::tables::{ANTILOG_TABLE, LOG_TABLE};

/// Multiplication in GF(256) using log tables.
pub const fn mul(a: u8, b: u8) -> u8 {
    if a == 0 || b == 0 {
        return 0;
    }

    let log_a = LOG_TABLE[a as usize] as u16;
    let log_b = LOG_TABLE[b as usize] as u16;
    let sum = (log_a + log_b) % 255;

    ANTILOG_TABLE[sum as usize].get()
}
