#![allow(clippy::module_name_repetitions)]

use std::num::TryFromIntError;

use ps_buffer::BufferError;
use thiserror::Error;

use crate::long::{LongEccHeaderConstructorError, LongEccHeaderFromByteSliceError};

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum GFError {
    #[error("Division by zero is undefined.")]
    DivByZero,
}

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialFromSliceError {
    #[error("Slice was {size} bytes, max 255 allowed in GF(256).")]
    TooLong { size: usize },
}

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialSetCoefficientsError {
    #[error("Range {offset}..{end} exceeds maximum coefficient index 254.")]
    OutOfBounds { offset: u8, end: usize },
}

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialMulError {
    #[error("Result degree exceeds maximum of 254.")]
    DegreeOverflow,
}

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialDivError {
    #[error("Division by zero polynomial.")]
    ZeroDivisor,
    #[error(transparent)]
    GFError(#[from] GFError),
}

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialXorError {
    #[error("Coefficient iterator exceeded maximum length of 255.")]
    TooManyCoefficients,
}

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum EuclideanError {
    /// The error-correction capability `t` requires the polynomial `x^(2t)`,
    /// whose degree exceeds [`Polynomial::MAX_DEGREE`](crate::Polynomial::MAX_DEGREE).
    #[error("Error-correction capability {t} exceeds the maximum of 127.")]
    CapabilityTooHigh { t: u8 },
    #[error(transparent)]
    PolynomialDiv(#[from] PolynomialDivError),
    #[error(transparent)]
    PolynomialMul(#[from] PolynomialMulError),
}

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum RSConstructorError {
    #[error("Parity count must be <= 63.")]
    ParityTooHigh,
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSGenerateParityError {
    #[error(transparent)]
    Division(#[from] PolynomialDivError),
    #[error(transparent)]
    SetCoefficients(#[from] PolynomialSetCoefficientsError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSEncodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    RSGenerateParityError(#[from] RSGenerateParityError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSComputeErrorsError {
    #[error(transparent)]
    GFError(#[from] GFError),
    #[error(transparent)]
    EuclideanError(#[from] EuclideanError),
    #[error("Too many errors; input is unrecoverable.")]
    TooManyErrors,
    #[error("The error locator derivative evaluated to zero.")]
    ZeroErrorLocatorDerivative,
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSDecodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error("Input length {received} is less than the parity length {parity_bytes}.")]
    InsufficientLength { parity_bytes: u8, received: usize },
    #[error(transparent)]
    RSComputeErrorsError(#[from] RSComputeErrorsError),
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    #[error("Too many errors to correct. Error computation nevertheless returned a valid polynomial, which is unlikely. Usually you'll get RSComputeErrorsError(TooManyErrors) instead.")]
    TooManyErrors,
    #[error(transparent)]
    TryFromIntError(#[from] TryFromIntError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum EncodeError {
    #[error(transparent)]
    LongEccEncodeError(#[from] LongEccEncodeError),
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    #[error(transparent)]
    RSEncodeError(#[from] RSEncodeError),
}

#[derive(Error, Debug, Clone)]
pub enum DecodeError {
    #[error("Insufficient input bytes for parity count of {0}: {0} * 2 > {1}.")]
    InsufficientParityBytes(u8, u8),
    #[error(transparent)]
    LongEccDecodeError(#[from] LongEccDecodeError),
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    #[error(transparent)]
    RSDecodeError(#[from] RSDecodeError),
}
#[derive(Error, Debug, Clone)]
pub enum EccError {
    #[error(transparent)]
    EncodeError(#[from] EncodeError),
    #[error(transparent)]
    DecodeError(#[from] DecodeError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccEncodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error("Parity {0} >= 64, which is too high.")]
    InvalidParity(u8),
    #[error("Invalid segment-to-parity ratio: {0} <= 2 * {1}")]
    InvalidSegmentParityRatio(u8, u8),
    #[error("Long ECC Header construction failed: {0}")]
    LongEccHeaderCtor(#[from] LongEccHeaderConstructorError),
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    #[error(transparent)]
    RSGenerateParityError(#[from] RSGenerateParityError),
    #[error(transparent)]
    TryFromIntError(#[from] TryFromIntError),
}

#[derive(Error, Debug, Clone)]
pub enum LongEccDecodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error("Integrity check failed after correction")]
    IntegrityCheckFailed,
    #[error("Codeword is invalid.")]
    InvalidCodeword,
    #[error("Failed to decode header: {0}")]
    HeaderDecode(#[from] LongEccHeaderFromByteSliceError),
    #[error("Failed to read data bytes.")]
    ReadDataError,
    #[error("Failed to read parity bytes.")]
    ReadParityError,
    #[error(transparent)]
    RSDecodeError(#[from] RSDecodeError),
    #[error(transparent)]
    TryFromIntError(#[from] TryFromIntError),
}
