use crate::long::LongEccHeader;

impl LongEccHeader {
    /// Serializes the header into its 32-byte wire format.
    #[must_use]
    pub fn to_bytes(self) -> [u8; 32] {
        let mut bytes = [0; 32];

        bytes[0..2].copy_from_slice(&self.magic.to_be_bytes());
        bytes[2] = self.version;
        bytes[3] = self.parity;
        bytes[4..8].copy_from_slice(&self.full_length.to_be_bytes());
        bytes[8..12].copy_from_slice(&self.message_length.to_be_bytes());
        bytes[12..16].copy_from_slice(&self.header_checksum.to_be_bytes());
        bytes[16..24].copy_from_slice(&self.checksum.to_be_bytes());
        bytes[24..32].copy_from_slice(&self.header_parity);

        bytes
    }
}

#[cfg(test)]
mod tests {
    use crate::long::{LongEccHeader, OverlapFactor, LONG_ECC_HEADER_MAGIC};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_to_bytes_layout() -> Result<(), TestError> {
        let checksum = 0x0123_4567_89AB_CDEF;
        let header = LongEccHeader::new(3, OverlapFactor::Double, 138, 100, checksum)?;

        let bytes = header.to_bytes();

        assert_eq!(bytes[0..2], LONG_ECC_HEADER_MAGIC.to_be_bytes());
        assert_eq!(bytes[2], 1);
        assert_eq!(bytes[3], 3 | OverlapFactor::Double.to_u8());
        assert_eq!(bytes[4..8], 138u32.to_be_bytes());
        assert_eq!(bytes[8..12], 100u32.to_be_bytes());
        assert_eq!(bytes[12..16], header.header_checksum().to_be_bytes());
        assert_eq!(bytes[16..24], checksum.to_be_bytes());
        assert_eq!(bytes[24..32], header.header_parity());

        Ok(())
    }

    #[test]
    fn test_to_bytes_from_bytes_roundtrip() -> Result<(), TestError> {
        let header = LongEccHeader::new(5, OverlapFactor::Quadruple, 482, 400, 42)?;

        let reparsed = LongEccHeader::from_bytes(header.to_bytes())?;

        assert_eq!(reparsed, header);

        Ok(())
    }
}
