use crate::{LongEccDecodeError, ReedSolomon, MAX_PARITY_BYTES};

use super::checksums::xxh64;
use super::fast_validate::fast_validate;
use super::{LongEccHeader, HEADER_SIZE};

/// Correct errors in-place in a codeword
pub fn correct_in_place(codeword: &mut [u8]) -> Result<LongEccHeader, LongEccDecodeError> {
    use LongEccDecodeError::{
        IntegrityCheckFailed, InvalidCodeword, ReadDataError, ReadParityError,
    };

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

    let excessive_parity = parity_bytes >= segment_distance.min(MAX_PARITY_BYTES as usize);
    let excessive_last_segment_length = last_segment_length > segment_length;

    if excessive_parity || excessive_last_segment_length {
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

    // Reed-Solomon can miscorrect a segment carrying more errors than its parity
    // can fix, landing on a wrong-but-valid codeword. The hash is the only guard,
    // so re-verify it after correction rather than trusting the corrected bytes.
    if header.parity > 0 && xxh64(&codeword[HEADER_SIZE..]) != header.xxh64 {
        return Err(IntegrityCheckFailed);
    }

    Ok(header)
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::LongEccDecodeError;

    use super::super::{decode, encode, HEADER_SIZE};
    use super::correct_in_place;

    type TestError = Box<dyn std::error::Error>;

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
    fn test_correct_in_place_rejects_hash_mismatch() -> Result<(), TestError> {
        let message = b"Integrity guard after correction".to_buffer()?;
        let parity = 2;
        let segment_length = 16;
        let segment_distance = 12;

        let mut encoded = encode(&message, parity, segment_length, segment_distance)?;

        // Drive a genuine miscorrection: a segment with fewer than `t` nonzero
        // symbols Reed-Solomon-decodes to the all-zeros codeword. Zeroing the body
        // and leaving a single nonzero symbol makes the corrector actively rewrite
        // it to all zeros and report success, yet the corrected body cannot match
        // the hash stored in the (untouched) header. The guard must reject it.
        encoded[HEADER_SIZE..].fill(0);
        encoded[HEADER_SIZE] = 1;

        let result = correct_in_place(&mut encoded);

        assert!(matches!(
            result,
            Err(LongEccDecodeError::IntegrityCheckFailed)
        ));

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
