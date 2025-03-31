use crate::{DecodeError, EncodeError, ReedSolomon};

/// Encodes a message by prepending an error-correcting code.
/// # Errors
/// - `RSConstructorError` is returned if `len(message) + 2 * parity` > `255`.
/// - `RSEncodeError` is returned if encoding fails for any reason.
pub fn encode(message: &[u8], parity: u8) -> Result<Vec<u8>, EncodeError> {
    let rs = ReedSolomon::new(parity)?;

    Ok(rs.encode(message)?)
}

/// Verifies the error-correcting code and returns the message.
/// # Errors
/// - `InputTooLarge` is returned if `len(received)` > 255 bytes.
/// - `InsufficientParityBytes` is returned if `parity > length / 2`.
/// - `RSDecodeError` is returned if decoding fails for any reason.
pub fn decode(received: &[u8], parity: u8) -> Result<Vec<u8>, DecodeError> {
    let length = received.len();
    let length: u8 = length
        .try_into()
        .or(Err(DecodeError::InputTooLarge(u32::try_from(length)?)))?;

    if parity > length >> 1 {
        return Err(DecodeError::InsufficientParityBytes(parity, length));
    }

    let rs = ReedSolomon::new(parity)?;

    Ok(rs.decode(received)?)
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

            assert_eq!(test_str, decoded);
        }

        Ok(())
    }
}
