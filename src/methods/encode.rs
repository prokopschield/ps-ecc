use ps_buffer::Buffer;

use crate::{long, EncodeError, ReedSolomon};

/// Encodes a message by adding an error-correcting code.
///
/// Messages longer than `255 - 2 * parity` bytes are encoded in the long ECC
/// format; all other messages become a single Reed-Solomon codeword.
/// # Errors
/// For messages that fit a single codeword:
/// - `RSConstructorError` is returned if `parity` exceeds [`crate::MAX_PARITY`].
/// - `RSEncodeError` is returned if encoding fails.
///
/// For longer messages:
/// - `LongEccEncodeError` is returned if long encoding fails.
pub fn encode(message: &[u8], parity: u8) -> Result<Buffer, EncodeError> {
    if message.len() + (usize::from(parity) << 1) > 0xff {
        let codeword = long::encode(message, parity, long::OverlapFactor::Simple)?;

        return Ok(codeword);
    }

    let rs = ReedSolomon::new(parity)?;

    Ok(rs.encode(message)?)
}
