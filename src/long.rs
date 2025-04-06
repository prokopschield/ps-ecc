use std::ops::Add;

use ps_buffer::Buffer;

use crate::{LongEccEncodeError, ReedSolomon};

const HEADER_SIZE: usize = std::mem::size_of::<LongEccHeader>();

#[derive(Clone, Copy, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[allow(clippy::module_name_repetitions)]
#[repr(C, align(16))]
pub struct LongEccHeader {
    pub full_length: u32,
    pub message_length: u32,
    pub segment_count: u32,
    pub parity: u8,
    pub segment_length: u8,
    pub segment_distance: u8,
    pub last_segment_length: u8,
}

impl LongEccHeader {
    #[inline]
    pub fn to_bytes(self) -> [u8; HEADER_SIZE] {
        let mut bytes = [0u8; HEADER_SIZE];

        bytes[0x0..0x4].copy_from_slice(&self.full_length.to_le_bytes());
        bytes[0x4..0x8].copy_from_slice(&self.message_length.to_le_bytes());
        bytes[0x8..0xC].copy_from_slice(&self.segment_count.to_le_bytes());
        bytes[12] = self.parity;
        bytes[13] = self.segment_length;
        bytes[14] = self.segment_distance;
        bytes[15] = self.last_segment_length;

        bytes
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
        .saturating_sub(segment_distance.saturating_sub(1))
        .div_ceil(new_bytes_per_segment)
        .saturating_add(1);
    let full_length = base_len + parity_bytes * segment_count;
    let processed_length = full_length - parity_bytes;
    let last_segment_length = if processed_length % segment_distance == 0 {
        segment_distance
    } else {
        processed_length % segment_distance
    };

    header.full_length = u32::try_from(full_length)?;
    header.last_segment_length = u8::try_from(last_segment_length)?;
    header.segment_count = u32::try_from(segment_count)?;

    let mut codeword = Buffer::with_capacity(full_length)?;

    codeword.extend_from_slice(&header.to_bytes())?;
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
            break;
        }
    }

    Ok(codeword)
}
