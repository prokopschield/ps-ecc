use crate::{LongEccConstructorError, ReedSolomon};

use super::{LongEccHeader, HEADER_SIZE};

impl LongEccHeader {
    /// Parse header from bytes with error correction
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, LongEccConstructorError> {
        if bytes.len() < HEADER_SIZE {
            return Err(LongEccConstructorError::InsufficientHeaderBytes(
                bytes.len().try_into()?,
            ));
        }

        // Copy data to stack array for in-place correction
        let mut data = [0u8; 24];

        data.copy_from_slice(&bytes[0..24]);

        // Correct errors in header data using Reed-Solomon (no heap allocation)
        ReedSolomon::correct_detached_data_in_place(&bytes[24..32], &mut data)?;

        let header = Self {
            full_length: u32::from_le_bytes(data[0..4].try_into()?),
            message_length: u32::from_le_bytes(data[4..8].try_into()?),
            parity: data[8],
            segment_length: data[9],
            segment_distance: data[10],
            last_segment_length: data[11],
            crc32: u32::from_le_bytes(data[12..16].try_into()?),
            xxh64: u64::from_le_bytes(data[16..24].try_into()?),
        };

        Ok(header)
    }
}

#[cfg(test)]
mod tests {
    use super::LongEccHeader;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_long_ecc_header_from_bytes() -> Result<(), TestError> {
        let bytes = [
            0x10, 0x00, 0x00, 0x00, // full length
            0x08, 0x00, 0x00, 0x00, // message length
            0x04, 0x0A, 0x05, 0x0B, // parity, seglen, segdist, lastseglen
            0x00, 0x00, 0x00, 0x00, // crc32
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // xxh64
            0xe0, 0x20, 0x7e, 0x4f, 0xd1, 0xc0, 0xbf, 0xae, // parity bytes
        ];

        let header = LongEccHeader::from_bytes(&bytes)?;

        assert_eq!(header.full_length, 16);
        assert_eq!(header.message_length, 8);
        assert_eq!(header.parity, 4);
        assert_eq!(header.segment_length, 10);
        assert_eq!(header.segment_distance, 5);
        assert_eq!(header.last_segment_length, 11);

        Ok(())
    }
}
