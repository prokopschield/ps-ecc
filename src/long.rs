use std::ops::Add;

use ps_buffer::Buffer;

use crate::{
    codeword::Codeword, LongEccConstructorError, LongEccDecodeError, LongEccEncodeError,
    LongEccToBytesError, ReedSolomon,
};

const HEADER_SIZE: usize = std::mem::size_of::<LongEccHeader>();

#[derive(Clone, Copy, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[allow(clippy::module_name_repetitions)]
#[repr(C, align(16))]
pub struct LongEccHeader {
    pub full_length: u32,
    pub message_length: u32,
    pub parity: u8,
    pub segment_length: u8,
    pub segment_distance: u8,
    pub last_segment_length: u8,
    pub crc32: u32, // CRC32 of message bytes
    pub xxh64: u64, // XXH64 of message + parity
}

impl LongEccHeader {
    /// Parse header from bytes with error correction
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, LongEccConstructorError> {
        if bytes.len() < HEADER_SIZE {
            return Err(LongEccConstructorError::InsufficientHeaderBytes(
                bytes.len().try_into()?,
            ));
        }

        // Extract data and parity portions
        let data_bytes = &bytes[0..24];
        let parity_bytes = &bytes[24..32];

        // Correct errors in header data using Reed-Solomon
        let corrected_data = ReedSolomon::correct_detached(parity_bytes, data_bytes)?;

        let header = Self {
            full_length: u32::from_le_bytes(corrected_data[0..4].try_into()?),
            message_length: u32::from_le_bytes(corrected_data[4..8].try_into()?),
            parity: corrected_data[8],
            segment_length: corrected_data[9],
            segment_distance: corrected_data[10],
            last_segment_length: corrected_data[11],
            crc32: u32::from_le_bytes(corrected_data[12..16].try_into()?),
            xxh64: u64::from_le_bytes(corrected_data[16..24].try_into()?),
        };

        Ok(header)
    }

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

/// Fast validation using checksums
pub fn fast_validate(codeword: &[u8]) -> Result<bool, LongEccDecodeError> {
    if codeword.len() < HEADER_SIZE {
        return Ok(false);
    }

    let header = LongEccHeader::from_bytes(codeword)?;

    // Validate CRC32 of message
    let message_start = HEADER_SIZE;
    let message_end = HEADER_SIZE + header.message_length as usize;

    if message_end > codeword.len() {
        return Ok(false);
    }

    let message = &codeword[message_start..message_end];
    let calculated_crc32 = crc32(message);

    if calculated_crc32 != header.crc32 {
        return Ok(false);
    }

    // Validate XXH64 of message + parity

    let calculated_xxh64 = xxh64(&codeword[HEADER_SIZE..]);

    if calculated_xxh64 != header.xxh64 {
        return Ok(false);
    }

    Ok(true)
}

/// Calculate CRC32 checksum
fn crc32(data: &[u8]) -> u32 {
    let mut hasher = crc32fast::Hasher::new();
    hasher.update(data);
    hasher.finalize()
}

/// Calculate XXH64 checksum
fn xxh64(data: &[u8]) -> u64 {
    xxhash_rust::xxh64::xxh64(data, 8_418_112_963_040_338_442)
}

/// Encode a message with long ECC protection
pub fn encode(
    message: &[u8],
    parity: u8,
    segment_length: u8,
    segment_distance: u8,
) -> Result<Buffer, LongEccEncodeError> {
    use LongEccEncodeError::{InvalidParity, InvalidSegmentParityRatio};

    // Validate parameters
    if parity >= 64 {
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
        codeword.extend_from_slice(&parity_data)?;

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

/// Correct errors in-place in a codeword
pub fn correct_in_place(codeword: &mut [u8]) -> Result<LongEccHeader, LongEccDecodeError> {
    use LongEccDecodeError::{InvalidCodeword, ReadDataError, ReadParityError};

    // Fast path - skip correction if data is valid
    if matches!(fast_validate(codeword), Ok(true)) {
        return Ok(LongEccHeader::from_bytes(codeword)?);
    }

    // Parse and correct header
    let header = LongEccHeader::from_bytes(codeword)?;

    // Extract parameters
    let parity_bytes = usize::from(header.parity) << 1;
    let last_segment_length = usize::from(header.last_segment_length);
    let segment_length = usize::from(header.segment_length);
    let segment_distance = usize::from(header.segment_distance);

    let mut parity_index = codeword.len().saturating_sub(parity_bytes);
    let mut data_index = parity_index.saturating_sub(last_segment_length);

    if parity_bytes >= segment_distance.min(127) || last_segment_length > segment_length {
        return Err(InvalidCodeword);
    }

    // Correct last chunk (data + parity)
    let (md, mp) = codeword[data_index..].split_at_mut(last_segment_length);
    ReedSolomon::correct_detached_in_place(mp, md)?;

    // Correct previous segments in reverse order
    while data_index > HEADER_SIZE {
        // Move indices to previous segment
        data_index = data_index.saturating_sub(segment_distance);
        parity_index = parity_index.saturating_sub(parity_bytes);

        // Define ranges for data and parity
        let data_range = data_index..data_index + segment_length;
        let parity_range = ..parity_bytes;

        // Split buffer at parity index
        let (data, parity) = codeword.split_at_mut(parity_index);

        // Get mutable references to data and parity with bounds checking
        let parity = parity.get_mut(parity_range).ok_or(ReadParityError)?;
        let data = data.get_mut(data_range).ok_or(ReadDataError)?;

        // Correct this segment
        ReedSolomon::correct_detached_in_place(parity, data)?;
    }

    Ok(header)
}

/// Decode a codeword to extract the original message
pub fn decode(codeword: &[u8]) -> Result<Codeword, LongEccDecodeError> {
    let mut buffer = Buffer::from_slice(codeword)?;
    let header = correct_in_place(&mut buffer)?;

    let message_length = usize::try_from(header.message_length)?;
    let codeword = Codeword {
        codeword: buffer.into(),
        range: HEADER_SIZE..HEADER_SIZE + message_length,
    };

    Ok(codeword)
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use super::*;

    #[derive(thiserror::Error, Debug)]
    enum TestError {
        #[error(transparent)]
        LongEccConstructor(#[from] LongEccConstructorError),
        #[error(transparent)]
        LongEccEncode(#[from] LongEccEncodeError),
        #[error(transparent)]
        LongEccDecode(#[from] LongEccDecodeError),
        #[error(transparent)]
        LongEccToBytes(#[from] LongEccToBytesError),
        #[error(transparent)]
        RSConstructorError(#[from] crate::RSConstructorError),
        #[error(transparent)]
        RSGenerateParityError(#[from] crate::RSGenerateParityError),
        #[error(transparent)]
        Buffer(#[from] ps_buffer::BufferError),
    }

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
    fn test_long_ecc_correct_in_place_no_errors() -> Result<(), TestError> {
        let message = b"Correct No Errors".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        let header = correct_in_place(&mut encoded)?;
        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );
        Ok(())
    }

    #[test]
    fn test_long_ecc_correct_in_place_one_error() -> Result<(), TestError> {
        let message = b"Correct One Error".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        encoded[HEADER_SIZE + 5] ^= 0b0000_0001;
        let header = correct_in_place(&mut encoded)?;
        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );
        Ok(())
    }

    #[test]
    fn test_long_ecc_correct_in_place_error_in_parity() -> Result<(), TestError> {
        let message = b"Error In Parity".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        let parity_start = HEADER_SIZE + message.len();
        encoded[parity_start + 1] ^= 0b0000_0010;
        let header = correct_in_place(&mut encoded)?;
        assert_eq!(header.message_length as usize, message.len());
        // We can't directly check the parity bytes, but if decode works, it's likely the parity was corrected.
        let decoded = decode(&encoded)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_long_ecc_correct_in_place_multiple_errors_recoverable() -> Result<(), TestError> {
        let message = b"Multiple Recoverable".to_buffer()?;
        let parity: u8 = 3;
        let segment_length: u8 = 12;
        let segment_distance: u8 = 8;
        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        encoded[HEADER_SIZE + 1] ^= 0b0000_0001;
        encoded[HEADER_SIZE + 7] ^= 0b0000_0010;
        encoded[HEADER_SIZE + message.len() + 3] ^= 0b0000_0100;
        let header = correct_in_place(&mut encoded)?;
        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );
        Ok(())
    }

    #[test]
    fn test_long_ecc_correct_in_place_too_many_errors() -> Result<(), TestError> {
        let message = b"Too Many Errors".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        encoded[HEADER_SIZE + 1] ^= 0b0000_0001;
        encoded[HEADER_SIZE + 3] ^= 0b0000_0010;
        encoded[HEADER_SIZE + 5] ^= 0b0000_0100; // More errors than can be corrected by parity=2
        let result = correct_in_place(&mut encoded);
        assert!(matches!(
            result,
            Err(LongEccDecodeError::RSDecodeError(
                crate::RSDecodeError::RSComputeErrorsError(
                    crate::RSComputeErrorsError::TooManyErrors
                )
            ))
        ));
        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_no_errors() -> Result<(), TestError> {
        let message = b"Decode No Errors".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let encoded = encode(&message, parity, segment_length, segment_distance)?;
        let decoded = decode(&encoded)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_one_error() -> Result<(), TestError> {
        let message = b"Decode One Error".to_buffer()?;
        let parity: u8 = 2;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        encoded[HEADER_SIZE + 2] ^= 0b0000_0001;
        let decoded = decode(&encoded)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_too_many_errors() -> Result<(), TestError> {
        let message = b"Decode Too Many".to_buffer()?;
        let parity: u8 = 1;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        encoded[HEADER_SIZE + 1] ^= 0b0000_0001;
        encoded[HEADER_SIZE + 3] ^= 0b0000_0010;
        let result = decode(&encoded);
        assert!(matches!(
            result,
            Err(LongEccDecodeError::RSDecodeError(
                crate::RSDecodeError::RSComputeErrorsError(
                    crate::RSComputeErrorsError::TooManyErrors
                )
            ))
        ));
        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_truncated_data() -> Result<(), TestError> {
        let header = LongEccHeader {
            full_length: 20,
            message_length: 10,
            parity: 1,
            segment_length: 8,
            segment_distance: 6,
            last_segment_length: 8,
            crc32: 0,
            xxh64: 0,
        };
        let mut truncated = header.to_bytes()?.to_vec();
        truncated.extend_from_slice(&[0u8; 5]); // Only 5 message bytes instead of 10
        let result = decode(&truncated);
        eprintln!("{result:?}");
        assert!(matches!(
            result,
            Err(crate::LongEccDecodeError::RSDecodeError(
                crate::RSDecodeError::RSComputeErrorsError(
                    crate::RSComputeErrorsError::TooManyErrors
                )
            ))
        ));
        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_incorrect_full_length() -> Result<(), TestError> {
        let message = b"Wrong Length".to_buffer()?;
        let parity: u8 = 1;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let encoded = encode(&message, parity, segment_length, segment_distance)?;
        let mut incorrect_length = encoded.to_vec();
        incorrect_length.pop();
        let result = decode(&incorrect_length);
        assert!(result.is_err());
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
    fn test_long_ecc_correct_in_place_zero_parity() -> Result<(), TestError> {
        let message = b"Zero Parity Correct".to_buffer()?;
        let parity: u8 = 0;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        let header = correct_in_place(&mut encoded)?;
        assert_eq!(header.parity, 0);
        assert_eq!(&encoded[HEADER_SIZE..], &message[..]);
        Ok(())
    }

    #[test]
    fn test_long_ecc_decode_zero_parity() -> Result<(), TestError> {
        let message = b"Zero Parity Decode".to_buffer()?;
        let parity: u8 = 0;
        let segment_length: u8 = 10;
        let segment_distance: u8 = 8;
        let encoded = encode(&message, parity, segment_length, segment_distance)?;
        let decoded = decode(&encoded)?;
        assert_eq!(&decoded[..], &message[..]);
        Ok(())
    }
}

#[cfg(test)]
mod encode_refactor_tests {
    use super::*;
    use ps_buffer::ToBuffer;

    #[derive(thiserror::Error, Debug)]
    enum TestError {
        #[error(transparent)]
        LongEccConstructor(#[from] LongEccConstructorError),
        #[error(transparent)]
        LongEccEncode(#[from] LongEccEncodeError),
        #[error(transparent)]
        LongEccDecode(#[from] LongEccDecodeError),
        #[error(transparent)]
        LongEccToBytes(#[from] LongEccToBytesError),
        #[error(transparent)]
        Buffer(#[from] ps_buffer::BufferError),
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
        let parity = 63; // Maximum valid parity value
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

#[cfg(test)]
mod correct_in_place_tests {
    use super::*;
    use ps_buffer::ToBuffer;

    #[derive(thiserror::Error, Debug)]
    enum TestError {
        #[error(transparent)]
        LongEccConstructor(#[from] LongEccConstructorError),
        #[error(transparent)]
        LongEccEncode(#[from] LongEccEncodeError),
        #[error(transparent)]
        LongEccDecode(#[from] LongEccDecodeError),
        #[error(transparent)]
        LongEccToBytes(#[from] LongEccToBytesError),
        #[error(transparent)]
        Buffer(#[from] ps_buffer::BufferError),
    }

    #[test]
    fn test_correct_in_place_single_segment() -> Result<(), TestError> {
        let message = b"Single segment test".to_buffer()?;
        let parity = 2;
        let segment_length = 30;
        let segment_distance = 25;

        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );
        Ok(())
    }

    #[test]
    fn test_correct_in_place_multiple_segments() -> Result<(), TestError> {
        let message =
            b"This is a longer message that will span multiple segments for testing".to_buffer()?;
        let parity = 3;
        let segment_length = 20;
        let segment_distance = 15;

        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );
        Ok(())
    }

    #[test]
    fn test_correct_in_place_error_correction_in_middle_segment() -> Result<(), TestError> {
        let message =
            b"Error correction in middle segment test with sufficient length".to_buffer()?;
        let parity = 2;
        let segment_length = 25;
        let segment_distance = 20;

        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;

        // Introduce an error in what should be a middle segment
        let error_position = HEADER_SIZE + 30;
        if error_position < encoded.len() {
            encoded[error_position] ^= 0b0000_0001;
        }

        let header = correct_in_place(&mut encoded)?;
        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );
        Ok(())
    }

    #[test]
    fn test_correct_in_place_edge_case_two_segments() -> Result<(), TestError> {
        let message = b"Two segment edge case".to_buffer()?;
        let parity = 1;
        let segment_length = 15;
        let segment_distance = 12;

        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;
        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );
        Ok(())
    }

    #[test]
    fn test_correct_in_place_with_parity_errors() -> Result<(), TestError> {
        let message = b"Parity error correction test".to_buffer()?;
        let parity = 2;
        let segment_length = 20;
        let segment_distance = 16;

        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;

        // Introduce errors in parity bytes of first segment
        let parity_start = HEADER_SIZE + message.len();
        if parity_start + 1 < encoded.len() {
            encoded[parity_start] ^= 0b0000_0001;
            encoded[parity_start + 1] ^= 0b0000_0010;
        }

        let header = correct_in_place(&mut encoded)?;
        assert_eq!(header.message_length as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );
        Ok(())
    }
}
