use crate::{DecodeError, EncodeError, ReedSolomon};

/// Encodes a message by prepending an error-correcting code.
/// # Errors
/// - `RSConstructorError` is returned if `len(message) + 2 * parity` > `255`.
/// - `RSEncodeError` is returned if encoding fails for any reason.
pub fn encode(message: &[u8], parity: usize) -> Result<Vec<u8>, EncodeError> {
    let k = message.len();
    let n = k + (parity << 1);
    let rs = ReedSolomon::new(n, k)?;

    Ok(rs.encode(message)?)
}

/// Verifies the error-correcting code and returns the message.
/// # Errors
/// - `RSConstructorError` is returned if `len(received)` > `255`.
/// - `RSConstructorError` is returned if `len(received)` > `255`.
/// - `RSEncodeError` is returned if encoding fails for any reason.
pub fn decode(received: &[u8], parity: usize) -> Result<Vec<u8>, DecodeError> {
    if parity > received.len() >> 1 {
        return Err(DecodeError::ParityTooLarge(parity, received.len()));
    }

    let n = received.len();
    let k = n - (parity << 1);
    let rs = ReedSolomon::new(n, k)?;

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
