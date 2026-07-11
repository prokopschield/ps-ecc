use std::num::NonZero;

use crate::error::GFError;

use super::tables::{ANTILOG_TABLE, LOG_TABLE};

/// Multiplicative inverse in GF(256).
pub const fn inv(a: u8) -> Result<NonZero<u8>, GFError> {
    use GFError::DivByZero;

    if a == 0 {
        return Err(DivByZero);
    }

    let value = ANTILOG_TABLE[(255 - LOG_TABLE[a as usize] as u16) as usize];

    Ok(value)
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use crate::GFError;

    #[test]
    fn result_size() {
        assert_eq!(std::mem::size_of::<Result<NonZero<u8>, GFError>>(), 1);
    }
}
