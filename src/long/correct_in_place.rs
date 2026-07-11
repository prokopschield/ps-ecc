use crate::{LongEccDecodeError, ReedSolomon};

use super::checksums::xxh64;
use super::fast_validate::fast_validate;
use super::{LongEccHeader, HEADER_SIZE};

/// Corrects errors in-place in a codeword.
pub fn correct_in_place(codeword: &mut [u8]) -> Result<LongEccHeader, LongEccDecodeError> {
    // Fast path - skip correction if data is valid
    if let Ok(Some(header)) = fast_validate(codeword) {
        return Ok(header);
    }

    correct_in_place_slow_path(codeword)
}

pub(crate) fn correct_in_place_slow_path(
    codeword: &mut [u8],
) -> Result<LongEccHeader, LongEccDecodeError> {
    use LongEccDecodeError::{
        IntegrityCheckFailed, InvalidCodeword, ReadDataError, ReadParityError,
    };

    // Parse and correct header
    let header = LongEccHeader::from_byte_slice(codeword)?;

    let full_length = usize::try_from(header.full_length())?;

    // Bytes beyond `full_length` are not part of the codeword; discard them
    // from correction and verification.
    let codeword = codeword.get_mut(..full_length).ok_or(InvalidCodeword)?;

    // Extract parameters
    let parity_bytes = usize::from(header.parity_bytes());
    let segment_length = usize::from(header.segment_length());
    let segment_distance = usize::from(header.segment_distance());

    // Reject headers whose parity leaves no room for new data within a segment;
    // the derived segment methods assume `2 * parity < segment_distance`.
    if parity_bytes >= segment_distance {
        return Err(InvalidCodeword);
    }

    let last_segment_length = usize::from(header.last_segment_length());

    // Reject codewords too short to contain the header, the final segment, and its parity.
    if codeword.len() < HEADER_SIZE + last_segment_length + parity_bytes {
        return Err(InvalidCodeword);
    }

    let mut parity_index = codeword.len() - parity_bytes;
    let mut data_index = parity_index - last_segment_length;

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
    // The checksum covers the parity as well, so miscorrected parity is also
    // caught, and with zero parity it is the sole corruption check.
    let payload = codeword
        .get(HEADER_SIZE..full_length)
        .ok_or(InvalidCodeword)?;

    if xxh64(payload) != header.checksum() {
        return Err(IntegrityCheckFailed);
    }

    // Restore the canonical header bytes, repairing any header corruption.
    codeword[..HEADER_SIZE].copy_from_slice(&header.to_bytes());

    Ok(header)
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::{LongEccDecodeError, MAX_PARITY};

    use super::super::{decode, encode, OverlapFactor, HEADER_SIZE};
    use super::correct_in_place;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_long_ecc_correct_in_place_no_errors() -> Result<(), TestError> {
        let message = b"Correct No Errors".to_buffer()?;
        let parity: u8 = 2;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());
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

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        encoded[HEADER_SIZE + 5] ^= 0b0000_0001;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());
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

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let parity_start = HEADER_SIZE + message.len();

        encoded[parity_start + 1] ^= 0b0000_0010;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());

        // We can't directly check the parity bytes, but if decode works, it's likely the parity was corrected.
        let decoded = decode(&encoded)?;

        assert_eq!(&decoded[..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_long_ecc_correct_in_place_multiple_errors_recoverable() -> Result<(), TestError> {
        let message = b"Multiple Recoverable".to_buffer()?;
        let parity: u8 = 3;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        encoded[HEADER_SIZE + 1] ^= 0b0000_0001;
        encoded[HEADER_SIZE + 7] ^= 0b0000_0010;
        encoded[HEADER_SIZE + message.len() + 3] ^= 0b0000_0100;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());
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

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

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

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.parity(), 0);
        assert_eq!(&encoded[HEADER_SIZE..], &message[..]);

        Ok(())
    }

    #[test]
    fn test_long_ecc_correct_in_place_zero_parity_detects_corruption() -> Result<(), TestError> {
        let message = b"Zero Parity Corruption".to_buffer()?;
        let parity: u8 = 0;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        encoded[HEADER_SIZE + 3] ^= 0b0000_0001;

        let result = correct_in_place(&mut encoded);

        assert!(matches!(
            result,
            Err(LongEccDecodeError::IntegrityCheckFailed)
        ));

        Ok(())
    }

    #[test]
    fn test_correct_in_place_single_segment() -> Result<(), TestError> {
        let message = b"Single segment test".to_buffer()?;
        let parity = 2;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );

        Ok(())
    }

    #[test]
    fn test_correct_in_place_multiple_segments() -> Result<(), TestError> {
        // With parity 3, a segment holds 249 data bytes; 710 bytes span three segments.
        let message = b"This is a longer message that will span multiple segments for testing"
            .repeat(10)
            .to_buffer()?;
        let parity = 3;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        encoded[HEADER_SIZE + 100] ^= 0b0000_0001;
        encoded[HEADER_SIZE + 400] ^= 0b0000_0010;
        encoded[HEADER_SIZE + 700] ^= 0b0000_0100;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );

        Ok(())
    }

    #[test]
    fn test_correct_in_place_error_correction_in_middle_segment() -> Result<(), TestError> {
        // With parity 2, a segment holds 251 data bytes; 504 bytes span three segments.
        let message = b"Error correction in middle segment test with sufficient length"
            .repeat(8)
            .to_buffer()?;
        let parity = 2;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        // Introduce an error in the middle segment
        encoded[HEADER_SIZE + 300] ^= 0b0000_0001;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );

        Ok(())
    }

    #[test]
    fn test_correct_in_place_edge_case_two_segments() -> Result<(), TestError> {
        // With parity 1, a segment holds 253 data bytes; 273 bytes span two segments.
        let message = b"Two segment edge case".repeat(13).to_buffer()?;
        let parity = 1;

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );

        Ok(())
    }

    #[test]
    fn test_correct_in_place_max_parity() -> Result<(), TestError> {
        let message = b"Maximum parity correction test".to_buffer()?;

        let mut encoded = encode(&message, MAX_PARITY, OverlapFactor::Simple)?;

        encoded[HEADER_SIZE + 10] ^= 0b0000_0001;
        encoded[HEADER_SIZE + 20] ^= 0b0000_0010;

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.parity(), MAX_PARITY);
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

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

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

        let mut encoded = encode(&message, parity, OverlapFactor::Simple)?;

        // Introduce errors in parity bytes of first segment
        let parity_start = HEADER_SIZE + message.len();

        if parity_start + 1 < encoded.len() {
            encoded[parity_start] ^= 0b0000_0001;
            encoded[parity_start + 1] ^= 0b0000_0010;
        }

        let header = correct_in_place(&mut encoded)?;

        assert_eq!(header.message_length() as usize, message.len());
        assert_eq!(
            &encoded[HEADER_SIZE..HEADER_SIZE + message.len()],
            &message[..]
        );

        Ok(())
    }
}
