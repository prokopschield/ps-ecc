use std::ops::Add;

use ps_buffer::Buffer;

use crate::{LongEccEncodeError, ReedSolomon, MAX_PARITY};

use super::checksums::{crc32, xxh64};
use super::{LongEccHeader, HEADER_SIZE};

/// Encode a message with long ECC protection
pub fn encode(
    message: &[u8],
    parity: u8,
    segment_length: u8,
    segment_distance: u8,
) -> Result<Buffer, LongEccEncodeError> {
    use LongEccEncodeError::{InvalidParity, InvalidSegmentParityRatio};

    // Validate parameters
    if parity > MAX_PARITY {
        return Err(InvalidParity(parity));
    }

    if parity >= (segment_distance >> 1) {
        return Err(InvalidSegmentParityRatio(segment_distance, parity));
    }

    // Ensure segment_length is at least as large as segment_distance
    let segment_distance_u8 = segment_distance;
    let segment_length_u8 = segment_length.max(segment_distance);

    // Calculate encoding parameters
    let base_len = message.len();
    let parity_bytes_per_segment = usize::from(parity << 1);
    let segment_distance = usize::from(segment_distance_u8);
    let segment_length = usize::from(segment_length_u8);

    // Calculate how many data bytes we can fit in each segment after parity
    let new_bytes_per_segment = segment_distance.saturating_sub(parity_bytes_per_segment);

    // Calculate number of segments needed
    let segment_count = base_len
        .saturating_sub(segment_length.saturating_sub(1))
        .div_ceil(new_bytes_per_segment)
        .saturating_add(1);

    // Calculate total size
    let full_length = HEADER_SIZE + base_len + parity_bytes_per_segment * segment_count;
    let processed_length = full_length - parity_bytes_per_segment - HEADER_SIZE;
    let n = (processed_length.saturating_sub(segment_length)).div_ceil(segment_distance);
    let last_segment_length = if processed_length >= n * segment_distance {
        processed_length - n * segment_distance
    } else {
        segment_length
    };

    // Calculate checksums
    let crc32_checksum = crc32(message);

    // Create header
    let header = LongEccHeader {
        full_length: u32::try_from(full_length)?,
        last_segment_length: u8::try_from(last_segment_length)?,
        message_length: message.len().try_into()?,
        parity,
        segment_length: segment_length_u8,
        segment_distance: segment_distance_u8,
        crc32: crc32_checksum,
        xxh64: 0, // Will be calculated after parity generation
    };

    // Initialize output buffer
    let mut codeword = Buffer::with_capacity(full_length)?;

    codeword.extend_from_slice(header.to_bytes()?)?;
    codeword.extend_from_slice(message)?;

    // Early return if no parity is requested
    if parity == 0 {
        return Ok(codeword);
    }

    // Generate parity for each segment
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
            debug_assert_eq!(
                header.last_segment_length as usize, segment_length_actual,
                "Segment length mismatch"
            );

            break;
        }
    }

    // Update XXH64 checksum
    let xxh64_checksum = xxh64(&codeword[HEADER_SIZE..]);

    // Rebuild header with correct XXH64
    let header = LongEccHeader {
        xxh64: xxh64_checksum,
        ..header
    };

    // Update header in codeword
    codeword[..HEADER_SIZE].copy_from_slice(&header.to_bytes()?);

    debug_assert_eq!(
        header.full_length as usize,
        codeword.len(),
        "Encoded length mismatch"
    );

    Ok(codeword)
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::{LongEccEncodeError, MAX_PARITY};

    use super::super::{decode, encode, LongEccHeader, HEADER_SIZE};

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_long_ecc_encode_no_parity() -> Result<(), TestError> {
        let message = b"No Parity".to_buffer()?;

        let encoded = encode(&message, 0, 10, 5)?;

        assert_eq!(encoded.len(), HEADER_SIZE + message.len());

        let header = LongEccHeader::from_bytes(&encoded)?;

        assert_eq!(header.parity, 0);
        assert_eq!(header.full_length as usize, encoded.len());
        assert_eq!(header.message_length as usize, message.len());

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_invalid_parity() -> Result<(), TestError> {
        let message = b"Invalid Parity".to_buffer()?;

        let result = encode(&message, 64, 10, 5);

        assert!(matches!(result, Err(LongEccEncodeError::InvalidParity(64))));

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_invalid_segment_parity_ratio() -> Result<(), TestError> {
        let message = b"Invalid Ratio".to_buffer()?;

        let result = encode(&message, 5, 10, 8);

        assert!(matches!(
            result,
            Err(LongEccEncodeError::InvalidSegmentParityRatio(8, 5))
        ));

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_with_parity() -> Result<(), TestError> {
        let message = b"With Parity".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        assert!(encoded.len() > HEADER_SIZE + message.len());

        let header = LongEccHeader::from_bytes(&encoded)?;

        assert_eq!(header.parity, parity);
        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(header.segment_length, segment_length);
        assert_eq!(header.segment_distance, segment_distance);
        assert_eq!(header.full_length as usize, encoded.len());

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_uneven_segments() -> Result<(), TestError> {
        let message = b"Uneven Segments".to_buffer()?;
        let parity: u8 = 1;
        let segment_length: u8 = 7;
        let segment_distance: u8 = 5;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let header = LongEccHeader::from_bytes(&encoded)?;

        assert_ne!(header.last_segment_length, segment_length);
        assert_eq!(header.full_length as usize, encoded.len());

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_zero_segment_distance() -> Result<(), TestError> {
        let message = b"Zero Distance".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 0;

        assert_eq!(
            encode(&message, parity, segment_length, segment_distance),
            Err(LongEccEncodeError::InvalidSegmentParityRatio(0, 2))
        );

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_segment_length_smaller_than_distance() -> Result<(), TestError> {
        let message = b"Small Segment".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 5;
        let segment_distance: u8 = 10;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let header = LongEccHeader::from_bytes(&encoded)?;

        assert_eq!(header.segment_length, segment_distance); // segment_length is max(segment_length, segment_distance)

        Ok(())
    }

    #[test]
    fn test_long_ecc_encode_empty_message() -> Result<(), TestError> {
        let message = b"".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        assert_eq!(encoded.len(), HEADER_SIZE + (usize::from(parity) << 1));

        let header = LongEccHeader::from_bytes(&encoded)?;

        assert_eq!(header.message_length, 0);
        assert_eq!(header.full_length as usize, encoded.len());

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_basic() -> Result<(), TestError> {
        let message = b"Hello, World!".to_buffer()?;
        let parity = 2;
        let segment_length = 10;
        let segment_distance = 8;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_no_parity() -> Result<(), TestError> {
        let message = b"No parity test".to_buffer()?;
        let parity = 0;
        let segment_length = 10;
        let segment_distance = 5;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_empty_message() -> Result<(), TestError> {
        let message = b"".to_buffer()?;
        let parity = 2;
        let segment_length = 10;
        let segment_distance = 8;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_single_byte() -> Result<(), TestError> {
        let message = b"X".to_buffer()?;
        let parity = 1;
        let segment_length = 8;
        let segment_distance = 4;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_large_message() -> Result<(), TestError> {
        let message = vec![0x42u8; 1000].to_buffer()?;
        let parity = 4;
        let segment_length = 50;
        let segment_distance = 40;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_segment_length_equals_distance() -> Result<(), TestError> {
        let message = b"Equal segments".to_buffer()?;
        let parity = 2;
        let segment_length = 10;
        let segment_distance = 10;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_segment_length_smaller_than_distance() -> Result<(), TestError>
    {
        let message = b"Small segments".to_buffer()?;
        let parity = 2;
        let segment_length = 5;
        let segment_distance = 10;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_zero_segment_distance() -> Result<(), TestError> {
        let message = b"Zero distance".to_buffer()?;
        let parity = 0;
        let segment_length = 10;
        let segment_distance = 0;

        // This should fail with InvalidSegmentParityRatio
        let result = encode(&message, parity, segment_length, segment_distance);

        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_max_parity() -> Result<(), TestError> {
        let message = b"Max parity".to_buffer()?;
        let parity = 32; // Below the 64 limit but still high
        let segment_length = 69;
        let segment_distance = 67;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_min_parity() -> Result<(), TestError> {
        let message = b"Min parity".to_buffer()?;
        let parity = 1;
        let segment_length = 10;
        let segment_distance = 8;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_exact_segment_fit() -> Result<(), TestError> {
        // Message size that fits exactly into segments
        let message = vec![0x42u8; 32].to_buffer()?; // 32 bytes message
        let parity = 2;
        let segment_length = 10;
        let segment_distance = 8;

        // With 8 bytes data per segment (10-2*2), 32 bytes needs exactly 4 segments

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_uneven_last_segment() -> Result<(), TestError> {
        let message = b"This message will create an uneven last segment".to_buffer()?;
        let parity = 3;
        let segment_length = 15;
        let segment_distance = 12;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_header_only_message() -> Result<(), TestError> {
        // Message that's smaller than header size
        let message = b"Hi".to_buffer()?;
        let parity = 1;
        let segment_length = 20;
        let segment_distance = 15;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_high_parity_ratio() -> Result<(), TestError> {
        let message = b"High parity ratio".to_buffer()?;
        let parity = 3;
        let segment_length = 10;
        let segment_distance = 7; // parity (3) >= segment_distance (7) >> 1 (3)

        // This should fail with InvalidSegmentParityRatio
        let result = encode(&message, parity, segment_length, segment_distance);

        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_encode_refactor_roundtrip_maximum_values() -> Result<(), TestError> {
        let message = vec![0xFFu8; 100].to_buffer()?;
        let parity = MAX_PARITY;
        let segment_length = 255; // Near maximum u8 value
        let segment_distance = 128;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_consistent_header_fields() -> Result<(), TestError> {
        let message = b"Header field consistency check".to_buffer()?;
        let parity = 2;
        let segment_length = 15;
        let segment_distance = 12;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let header = LongEccHeader::from_bytes(&encoded)?;

        // Check header fields are correctly set
        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(header.parity, parity);
        assert_eq!(header.segment_length, segment_length.max(segment_distance));
        assert_eq!(header.segment_distance, segment_distance);
        assert_eq!(header.full_length as usize, encoded.len());

        // Decode and verify message integrity
        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_segment_count_calculation() -> Result<(), TestError> {
        // Test case where the original and refactored versions might differ in segment count calculation
        let message = vec![0xAAu8; 50].to_buffer()?;
        let parity = 2;
        let segment_length = 12;
        let segment_distance = 10;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_encode_refactor_edge_case_small_segment_distance() -> Result<(), TestError> {
        let message = b"Small segment distance test".to_buffer()?;
        let parity = 1;
        let segment_length = 50;
        let segment_distance = 4;

        let encoded = encode(&message, parity, segment_length, segment_distance)?;

        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }
}
