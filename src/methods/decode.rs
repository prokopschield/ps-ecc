use crate::{codeword::Codeword, long, DecodeError, ReedSolomon};

/// Verifies the error-correcting code and returns the message.
///
/// Correctable corruption is repaired. For codewords longer than 255 bytes,
/// the `parity` argument is ignored, since the header records the parity, and
/// bytes beyond the full length recorded in the header are discarded, so
/// input rejected by [`validate`](crate::validate) may still decode
/// successfully.
/// # Errors
/// For codewords of at most 255 bytes:
/// - `InsufficientParityBytes` is returned if `parity > length / 2`.
/// - `RSConstructorError` is returned if `parity` otherwise exceeds [`crate::MAX_PARITY`].
/// - `RSDecodeError` is returned if decoding fails.
///
/// For longer codewords:
/// - `LongEccDecodeError` is returned if decoding fails.
pub fn decode(received: &[u8], parity: u8) -> Result<Codeword<'_>, DecodeError> {
    if let Ok(length) = u8::try_from(received.len()) {
        if parity > length >> 1 {
            return Err(DecodeError::InsufficientParityBytes(parity, length));
        }

        let rs = ReedSolomon::new(parity)?;

        Ok(rs.decode(received)?)
    } else {
        Ok(long::decode(received)?)
    }
}

#[cfg(test)]
mod tests {
    use crate::{decode, encode, EccError};

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
