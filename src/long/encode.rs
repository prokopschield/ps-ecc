use std::ops::Add;

use ps_buffer::Buffer;

use crate::{LongEccEncodeError, ReedSolomon, MAX_PARITY};

use super::checksums::xxh64;
use super::{LongEccHeader, OverlapFactor, HEADER_SIZE};

/// Encode a message with long ECC protection
pub fn encode(
    message: &[u8],
    parity: u8,
    overlap_factor: OverlapFactor,
) -> Result<Buffer, LongEccEncodeError> {
    use LongEccEncodeError::{InvalidParity, InvalidSegmentParityRatio};

    // Validate parameters
    if parity > MAX_PARITY {
        return Err(InvalidParity(parity));
    }

    // Derive the segment geometry the header encodes: each segment is a full
    // RS(255) codeword, and consecutive segments start `segment_distance` apart.
    let parity_bytes_u8 = parity << 1;
    let segment_length_u8 = 255 - parity_bytes_u8;
    let segment_distance_u8 = segment_length_u8 / overlap_factor.count();

    if parity_bytes_u8 >= segment_distance_u8 {
        return Err(InvalidSegmentParityRatio(segment_distance_u8, parity));
    }

    // Calculate encoding parameters
    let base_len = message.len();
    let parity_bytes_per_segment = usize::from(parity_bytes_u8);
    let segment_distance = usize::from(segment_distance_u8);
    let segment_length = usize::from(segment_length_u8);

    // Calculate how many data bytes we can fit in each segment after parity
    let new_bytes_per_segment = segment_distance - parity_bytes_per_segment;

    // Calculate number of segments needed
    let segment_count = base_len
        .saturating_sub(segment_length - 1)
        .div_ceil(new_bytes_per_segment)
        .saturating_add(1);

    // Calculate total size
    let full_length = HEADER_SIZE + base_len + parity_bytes_per_segment * segment_count;

    // Initialize the output buffer, reserving zeroed space for the header;
    // the header is written last so that its checksum can cover the parity.
    let mut codeword = Buffer::with_capacity(full_length)?;

    codeword.extend_from_slice([0u8; HEADER_SIZE])?;
    codeword.extend_from_slice(message)?;

    let mut last_segment_length = 0;

    // Generate parity for each segment
    if parity > 0 {
        let rs = ReedSolomon::new(parity)?;

        let mut index: usize = HEADER_SIZE;

        loop {
            let segment_end = index.add(segment_length).min(codeword.len());
            let segment_length_actual = segment_end - index;

            // Generate and append parity for this segment
            let parity_data = rs.generate_parity(&codeword[index..segment_end])?;

            codeword.extend_from_slice(parity_data)?;

            // Move to next segment
            index += segment_distance;

            // Break if we've processed the last segment
            if segment_length_actual != segment_length {
                last_segment_length = segment_length_actual;

                break;
            }
        }
    }

    // Calculate the checksum of the message and parity
    let checksum = xxh64(&codeword[HEADER_SIZE..]);

    // Create header
    let header = LongEccHeader::new(
        parity,
        overlap_factor,
        u32::try_from(full_length)?,
        message.len().try_into()?,
        checksum,
    )?;

    codeword[..HEADER_SIZE].copy_from_slice(&header.to_bytes());

    debug_assert_eq!(header.segment_length() as usize, segment_length);
    debug_assert_eq!(header.segment_distance() as usize, segment_distance);
    debug_assert_eq!(header.segment_count() as usize, segment_count);

    debug_assert!(
        parity == 0 || header.last_segment_length() as usize == last_segment_length,
        "Segment length mismatch"
    );

    debug_assert_eq!(
        header.full_length() as usize,
        codeword.len(),
        "Encoded length mismatch"
    );

    Ok(codeword)
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::{LongEccEncodeError, MAX_PARITY};

    use super::super::{
        decode, encode, LongEccHeader, OverlapFactor, HEADER_SIZE, LONG_ECC_HEADER_MAGIC,
    };

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_long_ecc_encode_no_parity() -> Result<(), TestError> {
        let message = b"No Parity".to_buffer()?;

        let encoded = encode(&message, 0, OverlapFactor::Simple)?;

        assert_eq!(encoded.len(), HEADER_SIZE + message.len());

        let header = LongEccHeader::from_byte_slice(&encoded)?;

        assert_eq!(header.parity(), 0);
        assert_eq!(header.full_length() as usize, encoded.len());
        assert_eq!(header.message_length() as usize, message.len());

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_invalid_parity() -> Result<(), TestError> {
        let message = b"Invalid Parity".to_buffer()?;

        let result = encode(&message, 64, OverlapFactor::Simple);

        assert!(matches!(result, Err(LongEccEncodeError::InvalidParity(64))));

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_invalid_segment_parity_ratio() -> Result<(), TestError> {
        let message = b"Invalid Ratio".to_buffer()?;

        // With parity 32 and quadruple overlap, each segment holds 64 parity bytes,
        // but consecutive segments start only (255 - 64) / 4 = 47 bytes apart.
        let result = encode(&message, 32, OverlapFactor::Quadruple);

        assert!(matches!(
            result,
            Err(LongEccEncodeError::InvalidSegmentParityRatio(47, 32))
        ));

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_with_parity() -> Result<(), TestError> {
        let message = b"With Parity".to_buffer()?;
        let parity: u8 = 2;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        assert!(encoded.len() > HEADER_SIZE + message.len());

        let header = LongEccHeader::from_byte_slice(&encoded)?;

        assert_eq!(header.parity(), parity);
        assert_eq!(header.overlap_factor(), OverlapFactor::Simple);
        assert_eq!(header.message_length() as usize, message.len());
        assert_eq!(header.segment_length(), 255 - (parity << 1));
        assert_eq!(header.segment_distance(), header.segment_length());
        assert_eq!(header.full_length() as usize, encoded.len());

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_uneven_segments() -> Result<(), TestError> {
        let message = b"Uneven Segments".to_buffer()?;
        let parity: u8 = 1;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let header = LongEccHeader::from_byte_slice(&encoded)?;

        assert_ne!(header.last_segment_length(), header.segment_length());
        assert_eq!(header.full_length() as usize, encoded.len());

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_empty_message() -> Result<(), TestError> {
        let message = b"".to_buffer()?;
        let parity: u8 = 2;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        assert_eq!(encoded.len(), HEADER_SIZE + (usize::from(parity) << 1));

        let header = LongEccHeader::from_byte_slice(&encoded)?;

        assert_eq!(header.message_length(), 0);
        assert_eq!(header.full_length() as usize, encoded.len());

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_basic() -> Result<(), TestError> {
        let message = b"Hello, World!".to_buffer()?;
        let parity = 2;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_no_parity() -> Result<(), TestError> {
        let message = b"No parity test".to_buffer()?;
        let parity = 0;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_empty_message() -> Result<(), TestError> {
        let message = b"".to_buffer()?;
        let parity = 2;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_single_byte() -> Result<(), TestError> {
        let message = b"X".to_buffer()?;
        let parity = 1;

        let encoded = encode(&message, parity, OverlapFactor::Double)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_large_message() -> Result<(), TestError> {
        let message = vec![0x42u8; 1000].to_buffer()?;
        let parity = 4;

        let encoded = encode(&message, parity, OverlapFactor::Double)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_triple_overlap() -> Result<(), TestError> {
        let message = b"Each byte is covered by three codewords".to_buffer()?;
        let parity = 2;

        let encoded = encode(&message, parity, OverlapFactor::Triple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_quadruple_overlap() -> Result<(), TestError> {
        let message = b"Each byte is covered by four codewords".to_buffer()?;
        let parity = 1;

        let encoded = encode(&message, parity, OverlapFactor::Quadruple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_max_parity() -> Result<(), TestError> {
        let message = b"Max parity".to_buffer()?;
        let parity = 32; // Below the 64 limit but still high

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_min_parity() -> Result<(), TestError> {
        let message = b"Min parity".to_buffer()?;
        let parity = 1;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_exact_segment_fit() -> Result<(), TestError> {
        // Message size that fits exactly into segments: with parity 2, a segment
        // holds 251 data bytes, and consecutive segments contribute 247 new bytes.
        let message = vec![0x42u8; 2 * 247].to_buffer()?;
        let parity = 2;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_uneven_last_segment() -> Result<(), TestError> {
        let message = b"This message will create an uneven last segment".to_buffer()?;
        let parity = 3;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_header_only_message() -> Result<(), TestError> {
        // Message that's smaller than header size
        let message = b"Hi".to_buffer()?;
        let parity = 1;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_high_parity_ratio() -> Result<(), TestError> {
        let message = b"High parity ratio".to_buffer()?;

        // Maximum parity with triple overlap: 126 parity bytes per segment, but
        // consecutive segments start only (255 - 126) / 3 = 43 bytes apart.
        let result = encode(&message, MAX_PARITY, OverlapFactor::Triple);

        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_maximum_values() -> Result<(), TestError> {
        let message = vec![0xFFu8; 100].to_buffer()?;
        let parity = MAX_PARITY;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_consistent_header_fields() -> Result<(), TestError> {
        let message = b"Header field consistency check".to_buffer()?;
        let parity = 2;

        let encoded = encode(&message, parity, OverlapFactor::Double)?;

        let header = LongEccHeader::from_byte_slice(&encoded)?;

        // Check header fields are correctly set
        assert_eq!(header.magic(), LONG_ECC_HEADER_MAGIC);
        assert_eq!(header.version(), 1);
        assert_eq!(header.message_length() as usize, message.len());
        assert_eq!(header.parity(), parity);
        assert_eq!(header.overlap_factor(), OverlapFactor::Double);
        assert_eq!(header.segment_length(), 255 - (parity << 1));
        assert_eq!(header.segment_distance(), header.segment_length() / 2);
        assert_eq!(header.full_length() as usize, encoded.len());

        // The serialized header roundtrips through error correction unchanged.
        let reparsed = LongEccHeader::from_bytes(header.to_bytes())?;

        assert_eq!(reparsed.header_checksum(), header.header_checksum());
        assert_eq!(reparsed.header_parity(), header.header_parity());
        assert_eq!(reparsed, header);

        // Decode and verify message integrity
        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_segment_count_calculation() -> Result<(), TestError> {
        // With parity 4, a segment holds 247 data bytes and each full segment
        // contributes 239 new bytes, so 600 bytes require three segments.
        let message = vec![0xAAu8; 600].to_buffer()?;
        let parity = 4;

        let encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let header = LongEccHeader::from_byte_slice(&encoded)?;

        assert_eq!(header.segment_count(), 3);
        assert_eq!(
            encoded.len(),
            HEADER_SIZE + message.len() + 3 * usize::from(header.parity_bytes())
        );

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_edge_case_small_segment_distance() -> Result<(), TestError> {
        // Quadruple overlap yields the smallest segment distance: (255 - 2) / 4 = 63.
        let message = b"Small segment distance test".to_buffer()?;
        let parity = 1;

        let encoded = encode(&message, parity, OverlapFactor::Quadruple)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }
}
