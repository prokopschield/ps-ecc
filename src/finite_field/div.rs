use crate::error::GFError;

use super::inv::inv;
use super::mul::mul;

/// Division in GF(256).
pub fn div(a: u8, b: u8) -> Result<u8, GFError> {
    use GFError::DivByZero;

    if b == 0 {
        return Err(DivByZero);
    }

    if a == 0 {
        return Ok(0);
    }

    Ok(mul(a, inv(b)?.get()))
}
