use ps_buffer::Buffer;

use crate::error::RSEncodeError;
use crate::ReedSolomon;

impl ReedSolomon {
    /// Encodes a message into a codeword.
    /// # Errors
    /// - [`ps_buffer::BufferError`] is returned if memory allocation fails.
    /// - [`RSEncodeError::RSGenerateParityError`] if parity generation fails.
    pub fn encode(&self, message: &[u8]) -> Result<Buffer, RSEncodeError> {
        let mut buffer = Buffer::with_capacity(message.len() + usize::from(self.parity_bytes()))?;

        buffer.extend_from_slice(self.generate_parity(message)?)?;
        buffer.extend_from_slice(message)?;

        Ok(buffer)
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::ReedSolomon;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_encode_empty_message() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"".to_buffer()?;
        let encoded = rs.encode(&message)?;

        assert_eq!(encoded.len(), 4); // 2 parity * 2 bytes
        assert!(rs.validate(&encoded).is_none());

        Ok(())
    }
}
