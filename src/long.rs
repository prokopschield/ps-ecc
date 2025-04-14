use std::ops::Add;

use ps_buffer::Buffer;

use crate::{
    codeword::Codeword, LongEccConstructorError, LongEccDecodeError, LongEccEncodeError,
    ReedSolomon,
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
    pub fn to_bytes(self) -> [u8; HEADER_SIZE] {
        let mut bytes = [0u8; HEADER_SIZE];

        bytes[0x0..0x4].copy_from_slice(&self.full_length.to_le_bytes());
        bytes[0x4..0x8].copy_from_slice(&self.message_length.to_le_bytes());
        bytes[0x8] = self.parity;
        bytes[0x9] = self.segment_length;
        bytes[0xA] = self.segment_distance;
        bytes[0xB] = self.last_segment_length;

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

    let mut codeword = Buffer::with_capacity(full_length)?;

    codeword.extend_from_slice(header.to_bytes())?;
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

pub fn correct_in_place(codeword: &mut [u8]) -> Result<LongEccHeader, LongEccDecodeError> {
    use LongEccDecodeError::{ReadDataError, ReadParityError};

    let header = LongEccHeader::from_bytes(codeword)?;

    let parity_bytes = usize::from(header.parity) << 1;
    let last_segment_length = usize::from(header.last_segment_length);
    let segment_length = usize::from(header.segment_length);
    let segment_distance = usize::from(header.segment_distance);

    let mut parity_index = codeword.len().saturating_sub(parity_bytes);
    let mut data_index = parity_index.saturating_sub(last_segment_length);

    // last chunk
    let (md, mp) = codeword[data_index..].split_at_mut(last_segment_length);
    ReedSolomon::correct_both_detached_in_place(mp, md)?;

    while data_index > 0 {
        data_index = data_index.saturating_sub(segment_distance);
        parity_index = parity_index.saturating_sub(parity_bytes);

        let data_range = data_index..data_index + segment_length;
        let parity_range = ..parity_bytes;

        let (data, parity) = codeword.split_at_mut(parity_index);

        let parity = parity.get_mut(parity_range).ok_or(ReadParityError)?;
        let data = data.get_mut(data_range).ok_or(ReadDataError)?;

        ReedSolomon::correct_both_detached_in_place(parity, data)?;
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
