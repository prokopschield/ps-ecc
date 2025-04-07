use ps_buffer::Buffer;

use crate::{codeword::Codeword, long, DecodeError, EncodeError, ReedSolomon};

/// Encodes a message by adding an error-correcting code.
/// # Errors
/// - `RSConstructorError` is returned if `len(message) + 2 * parity` > `255`.
/// - `RSEncodeError` is returned if encoding fails for any reason.
pub fn encode(message: &[u8], parity: u8) -> Result<Buffer, EncodeError> {
    if message.len() + (usize::from(parity) << 1) > 0xff {
        let segment_length = 0xFF - (parity << 1);
        let codeword = long::encode(message, parity, segment_length, segment_length)?;

        return Ok(codeword);
    }

    let rs = ReedSolomon::new(parity)?;

    Ok(rs.encode(message)?)
}

/// Verifies the error-correcting code and returns the message.
/// # Errors
/// - `InputTooLarge` is returned if `len(received)` > 255 bytes.
/// - `InsufficientParityBytes` is returned if `parity > length / 2`.
/// - `RSDecodeError` is returned if decoding fails for any reason.
pub fn decode(received: &[u8], parity: u8) -> Result<Codeword, DecodeError> {
    if let Ok(length) = u8::try_from(received.len()) {
        if parity > length >> 1 {
            return Err(DecodeError::InsufficientParityBytes(parity, length));
        }

        let rs = ReedSolomon::new(parity)?;

        Ok(rs.decode(received)?.into())
    } else {
        Ok(long::decode(received)?)
    }
}

#[cfg(test)]
mod tests {
    use crate::EccError;

    use super::{decode, encode};

    #[test]
    fn ecc_works() -> Result<(), EccError> {
        let test_str = "Strč prst skrz krk! ¯\\_(ツ)_/¯".as_bytes();
        let mut encoded = encode(test_str, 13)?;

        for i in 0..13 {
            let index = (i * 37) % encoded.len();
            encoded[index] ^= (i * index + 13).to_le_bytes()[0];
            let decoded = decode(&encoded, 13)?;

            assert_eq!(test_str, &decoded[..]);
        }

        Ok(())
    }
}
