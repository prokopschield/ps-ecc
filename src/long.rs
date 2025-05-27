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
}

impl LongEccHeader {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, LongEccConstructorError> {
        if bytes.len() < HEADER_SIZE {
            return Err(LongEccConstructorError::InsufficientHeaderBytes(
                bytes.len().try_into()?,
            ));
        }

        let bytes = ReedSolomon::correct_detached(&bytes[12..16], &bytes[0..12])?;

        let header = Self {
            full_length: u32::from_le_bytes(bytes[0..4].try_into()?),
            message_length: u32::from_le_bytes(bytes[4..8].try_into()?),
            parity: *bytes.get(8).unwrap_or(&0),
            segment_length: *bytes.get(9).unwrap_or(&0),
            segment_distance: *bytes.get(10).unwrap_or(&0),
            last_segment_length: *bytes.get(11).unwrap_or(&0),
        };

        Ok(header)
    }

    #[inline]
    pub fn to_bytes(self) -> Result<[u8; HEADER_SIZE], LongEccToBytesError> {
        let mut bytes = [0u8; HEADER_SIZE];

        bytes[0x0..0x4].copy_from_slice(&self.full_length.to_le_bytes());
        bytes[0x4..0x8].copy_from_slice(&self.message_length.to_le_bytes());
        bytes[0x8] = self.parity;
        bytes[0x9] = self.segment_length;
        bytes[0xA] = self.segment_distance;
        bytes[0xB] = self.last_segment_length;

        let parity = ReedSolomon::new(2)?.generate_parity(&bytes[0..12])?;

        bytes[0xC..=0xF].copy_from_slice(&parity);

        Ok(bytes)
    }
}

pub fn encode(
    message: &[u8],
    parity: u8,
    segment_length: u8,
    segment_distance: u8,
) -> Result<Buffer, LongEccEncodeError> {
    use LongEccEncodeError::{InvalidParity, InvalidSegmentParityRatio};

    if parity >= 64 {
        return Err(InvalidParity(parity));
    }

    if parity >= (segment_distance >> 1) {
        return Err(InvalidSegmentParityRatio(segment_distance, parity));
    }

    let segment_length = segment_length.max(segment_distance);

    let mut header = LongEccHeader {
        message_length: message.len().try_into()?,
        parity,
        segment_length,
        segment_distance,
        ..Default::default()
    };

    let base_len = HEADER_SIZE + message.len();
    let parity_bytes = usize::from(parity << 1);
    let segment_distance = usize::from(segment_distance);
    let segment_length = usize::from(segment_length);
    let new_bytes_per_segment = segment_distance - parity_bytes;
    let segment_count = base_len
        .saturating_sub(segment_length.saturating_sub(1))
        .div_ceil(new_bytes_per_segment)
        .saturating_add(1);
    let full_length = base_len + parity_bytes * segment_count;
    let processed_length = full_length - parity_bytes;
    let n = (processed_length - segment_length).div_ceil(segment_distance);
    let last_segment_length = if processed_length >= n * segment_distance {
        processed_length - n * segment_distance
    } else {
        segment_length
    };

    header.full_length = u32::try_from(full_length)?;
    header.last_segment_length = u8::try_from(last_segment_length)?;

    let mut codeword = Buffer::with_capacity(full_length)?;

    codeword.extend_from_slice(header.to_bytes()?)?;
    codeword.extend_from_slice(message)?;

    if parity == 0 {
        return Ok(codeword);
    }

    let rs = ReedSolomon::new(parity)?;

    let mut index: usize = 0;

    loop {
        let next_segment = index..index.add(segment_length).min(codeword.len());
        let next_segment_length = next_segment.end - next_segment.start;

        codeword.extend_from_slice(&rs.generate_parity(&codeword[next_segment])?)?;
        index += segment_distance;

        if next_segment_length != segment_length {
            debug_assert_eq!(
                header.last_segment_length as usize, next_segment_length,
                "Next segment length doesn't match"
            );

            break;
        }
    }

    debug_assert_eq!(
        header.full_length as usize,
        codeword.len(),
        "Full length doesn't match"
    );

    Ok(codeword)
}

pub fn correct_in_place(codeword: &mut [u8]) -> Result<LongEccHeader, LongEccDecodeError> {
    use LongEccDecodeError::{InvalidCodeword, ReadDataError, ReadParityError};

    let header = LongEccHeader::from_bytes(codeword)?;

    let parity_bytes = usize::from(header.parity) << 1;
    let last_segment_length = usize::from(header.last_segment_length);
    let segment_length = usize::from(header.segment_length);
    let segment_distance = usize::from(header.segment_distance);

    let mut parity_index = codeword.len().saturating_sub(parity_bytes);
    let mut data_index = parity_index.saturating_sub(last_segment_length);

    if parity_bytes >= segment_distance.min(127) || last_segment_length > segment_length {
        return Err(InvalidCodeword);
    }

    // last chunk
    let (md, mp) = codeword[data_index..].split_at_mut(last_segment_length);
    ReedSolomon::correct_detached_in_place(mp, md)?;

    while data_index > 0 {
        data_index = data_index.saturating_sub(segment_distance);
        parity_index = parity_index.saturating_sub(parity_bytes);

        let data_range = data_index..data_index + segment_length;
        let parity_range = ..parity_bytes;

        let (data, parity) = codeword.split_at_mut(parity_index);

        let parity = parity.get_mut(parity_range).ok_or(ReadParityError)?;
        let data = data.get_mut(data_range).ok_or(ReadDataError)?;

        ReedSolomon::correct_detached_in_place(parity, data)?;
    }

    Ok(header)
}

pub fn decode(codeword: &[u8]) -> Result<Codeword, LongEccDecodeError> {
    let mut buffer = Buffer::from_slice(codeword)?;
    let header = correct_in_place(&mut buffer)?;
    let codeword = Codeword {
        codeword: buffer.into(),
        range: HEADER_SIZE..HEADER_SIZE + usize::try_from(header.message_length)?,
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
            0x2D, 0xE6, 0x64, 0x1A, // parity bytes
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
        };
        let bytes = header.to_bytes()?;
        assert_eq!(&bytes[0..4], &32u32.to_le_bytes());
        assert_eq!(&bytes[4..8], &12u32.to_le_bytes());
        assert_eq!(bytes[8], 6);
        assert_eq!(bytes[9], 15);
        assert_eq!(bytes[10], 7);
        assert_eq!(bytes[11], 18);
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
            )) // HEADER_SIZE + message_length = 12 + 10 = 22? No, 12 + 10 = 22. Header is 12. 12 + 10 = 22. Got 17.
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
        assert_eq!(encoded.len(), HEADER_SIZE + 3 * (usize::from(parity) << 1));
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
