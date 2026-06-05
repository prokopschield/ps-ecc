use crate::{LongEccToBytesError, ReedSolomon};

use super::{LongEccHeader, HEADER_SIZE};

impl LongEccHeader {
    /// Serialize header to bytes with Reed-Solomon parity
    #[inline]
    pub fn to_bytes(self) -> Result<[u8; HEADER_SIZE], LongEccToBytesError> {
        let mut bytes = [0u8; HEADER_SIZE];

        // Serialize header fields
        bytes[0..4].copy_from_slice(&self.full_length.to_le_bytes());
        bytes[4..8].copy_from_slice(&self.message_length.to_le_bytes());
        bytes[8] = self.parity;
        bytes[9] = self.segment_length;
        bytes[10] = self.segment_distance;
        bytes[11] = self.last_segment_length;
        bytes[12..16].copy_from_slice(&self.crc32.to_le_bytes());
        bytes[16..24].copy_from_slice(&self.xxh64.to_le_bytes());

        // Generate parity for header data
        let rs = ReedSolomon::new(4)?; // RS(32,8) - 24 data bytes + 8 parity bytes
        let parity = rs.generate_parity(&bytes[0..24])?;

        // Append parity bytes
        bytes[24..32].copy_from_slice(&parity);

        Ok(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::LongEccHeader;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_long_ecc_header_to_bytes() -> Result<(), TestError> {
        let header = LongEccHeader {
            full_length: 32,
            message_length: 12,
            parity: 6,
            segment_length: 15,
            segment_distance: 7,
            last_segment_length: 18,
            crc32: 0x1234_5678,
            xxh64: 0x1234_5678_90AB_CDEF,
        };

        let bytes = header.to_bytes()?;

        assert_eq!(&bytes[0..4], &32u32.to_le_bytes());
        assert_eq!(&bytes[4..8], &12u32.to_le_bytes());
        assert_eq!(bytes[8], 6);
        assert_eq!(bytes[9], 15);
        assert_eq!(bytes[10], 7);
        assert_eq!(bytes[11], 18);
        assert_eq!(&bytes[12..16], &0x1234_5678u32.to_le_bytes());
        assert_eq!(&bytes[16..24], &0x1234_5678_90AB_CDEFu64.to_le_bytes());

        Ok(())
    }
}
