use std::ops::Rem;

use crate::error::RSGenerateParityError;
use crate::{ParityBytes, Polynomial, ReedSolomon};

use super::super::generator::generator_poly;

impl ReedSolomon {
    /// Generates parity bytes.
    /// # Errors
    /// - [`RSGenerateParityError::SetCoefficients`] is returned if
    ///   `2 * parity + message.len()` exceeds 255.
    /// - [`RSGenerateParityError::Division`] is returned if the generator
    ///   polynomial is zero (not expected to occur).
    pub fn generate_parity(&self, message: &[u8]) -> Result<ParityBytes, RSGenerateParityError> {
        let mut p = Polynomial::default();

        p.set_coefficients(self.parity_bytes(), message)?;

        let remainder = p.rem(generator_poly(self.parity()))?;

        Ok(ParityBytes::new(&remainder, self.parity()))
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::ToBuffer;

    use crate::ReedSolomon;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_generate_parity_no_errors() -> Result<(), TestError> {
        let rs = ReedSolomon::new(4)?;
        let message = b"Test".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert_eq!(parity.len(), 8); // 4 parity * 2 bytes

        Ok(())
    }

    #[test]
    fn test_generate_parity_empty_message() -> Result<(), TestError> {
        let rs = ReedSolomon::new(2)?;
        let message = b"".to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert_eq!(parity.len(), 4);
        assert_eq!(parity.as_slice(), &[0, 0, 0, 0]);

        Ok(())
    }

    #[test]
    fn test_generate_parity_large_message() -> Result<(), TestError> {
        let rs = ReedSolomon::new(8)?;
        let message = vec![42; 200].to_buffer()?;
        let parity = rs.generate_parity(&message)?;

        assert_eq!(parity.len(), 16); // 8 parity * 2 bytes

        let encoded = rs.encode(&message)?;

        assert_eq!(&encoded[..16], parity.as_slice());

        Ok(())
    }
}
