use crate::long::LongEccHeader;

use super::LongEccHeaderFromBytesError;

#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccHeaderFromByteSliceError {
    #[error("Insufficient bytes for header: got {0}, need 32.")]
    InsufficientBytes(usize),

    #[error(transparent)]
    FromBytes(#[from] LongEccHeaderFromBytesError),
}

impl LongEccHeader {
    /// Parses a header from the first 32 bytes of a slice.
    ///
    /// # Errors
    /// Returns `InsufficientBytes` if `bytes` holds fewer than 32 bytes, or a
    /// [`LongEccHeaderFromBytesError`] if parsing the header fails.
    pub fn from_byte_slice(bytes: &[u8]) -> Result<Self, LongEccHeaderFromByteSliceError> {
        let len = bytes.len();

        let header_bytes = bytes
            .first_chunk::<32>()
            .ok_or(LongEccHeaderFromByteSliceError::InsufficientBytes(len))?;

        Ok(Self::from_bytes(*header_bytes)?)
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor};

    use super::LongEccHeaderFromByteSliceError;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_from_byte_slice_roundtrip() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Simple, 86, 50, 7)?;

        let parsed = LongEccHeader::from_byte_slice(&header.to_bytes())?;

        assert_eq!(parsed, header);

        Ok(())
    }

    #[test]
    fn test_from_byte_slice_ignores_trailing_bytes() -> Result<(), TestError> {
        let header = LongEccHeader::new(2, OverlapFactor::Double, 86, 50, 7)?;

        let mut bytes = header.to_bytes().to_vec();

        bytes.extend_from_slice(&[0xAB; 16]);

        let parsed = LongEccHeader::from_byte_slice(&bytes)?;

        assert_eq!(parsed, header);

        Ok(())
    }

    #[test]
    fn test_from_byte_slice_insufficient_bytes() {
        let result = LongEccHeader::from_byte_slice(&[0u8; 31]);

        assert_eq!(
            result,
            Err(LongEccHeaderFromByteSliceError::InsufficientBytes(31))
        );
    }

    #[test]
    fn test_from_byte_slice_empty_slice() {
        let result = LongEccHeader::from_byte_slice(&[]);

        assert_eq!(
            result,
            Err(LongEccHeaderFromByteSliceError::InsufficientBytes(0))
        );
    }
}
