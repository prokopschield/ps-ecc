use crate::long::checksums::xxh64;

pub(in super::super) fn calculate_header_checksum(
    version: u8,
    parity: u8,
    full_length: u32,
    message_length: u32,
) -> u32 {
    let mut bytes = [0u8; 10];

    bytes[0] = version;
    bytes[1] = parity;
    bytes[2..6].copy_from_slice(&full_length.to_be_bytes());
    bytes[6..10].copy_from_slice(&message_length.to_be_bytes());

    let hash = xxh64(&bytes);

    // Fold the 64-bit hash into 32 bits, preserving the entropy of both halves.
    #[allow(clippy::cast_possible_truncation)]
    let folded = ((hash >> 32) ^ hash) as u32;

    folded
}

#[cfg(test)]
mod tests {
    use super::calculate_header_checksum;

    #[test]
    fn test_header_checksum_deterministic() {
        assert_eq!(
            calculate_header_checksum(1, 2, 100, 50),
            calculate_header_checksum(1, 2, 100, 50)
        );
    }

    #[test]
    fn test_header_checksum_sensitive_to_each_field() {
        let base = calculate_header_checksum(1, 2, 100, 50);

        assert_ne!(base, calculate_header_checksum(2, 2, 100, 50));
        assert_ne!(base, calculate_header_checksum(1, 3, 100, 50));
        assert_ne!(base, calculate_header_checksum(1, 2, 101, 50));
        assert_ne!(base, calculate_header_checksum(1, 2, 100, 51));
    }

    #[test]
    fn test_header_checksum_nonzero_for_zero_fields() {
        // The former multiplicative combinator collapsed to 0 whenever any
        // factor was 0; the hash must not.
        assert_ne!(calculate_header_checksum(1, 0, 100, 50), 0);
        assert_ne!(calculate_header_checksum(1, 2, 100, 0), 0);
        assert_ne!(calculate_header_checksum(1, 0, 32, 0), 0);
    }

    #[test]
    fn test_header_checksum_distinguishes_swapped_lengths() {
        assert_ne!(
            calculate_header_checksum(1, 2, 100, 50),
            calculate_header_checksum(1, 2, 50, 100)
        );
    }
}
