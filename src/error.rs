#![allow(clippy::module_name_repetitions)]

use std::{array::TryFromSliceError, num::TryFromIntError};

use ps_buffer::BufferError;
use thiserror::Error;

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum GFError {
    #[error("Division by zero is undefined.")]
    DivByZero,
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum PolynomialError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    GFError(#[from] GFError),
    #[error("Divisor cannot be zero.")]
    ZeroDivisor,
}

#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum RSConstructorError {
    #[error("Parity count must be <= 127.")]
    ParityTooHigh,
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSGenerateParityError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    PolynomialError(#[from] PolynomialError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSEncodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    RSGenerateParityError(#[from] RSGenerateParityError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSValidationError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSEuclideanError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    PolynomialError(#[from] PolynomialError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSComputeErrorsError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    GFError(#[from] GFError),
    #[error(transparent)]
    RSEuclideanError(#[from] RSEuclideanError),
    #[error("Too many errors, input in unrecoverable.")]
    TooManyErrors,
    #[error("The error locator derivative evaluated to zero.")]
    ZeroErrorLocatorDerivative,
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSDecodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    RSComputeErrorsError(#[from] RSComputeErrorsError),
    #[error(transparent)]
    RSValidationError(#[from] RSValidationError),
    #[error("Too many errors to correct.")]
    TooManyErrors,
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

#[derive(Error, Debug, Clone, Copy)]
pub enum LongEccConstructorError {
    #[error(transparent)]
    TryFromSliceError(#[from] TryFromSliceError),
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccEncodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error("Parity {0} >= 64, which is too high.")]
    InvalidParity(u8),
    #[error("Invalid segment-to-parity ratio: {0} < 2 * {1}")]
    InvalidSegmentParityRatio(u8, u8),
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
    #[error("Codeword is invalid.")]
    InvalidCodeword,
    #[error(transparent)]
    LongEccConstructorError(#[from] LongEccConstructorError),
    #[error("Failed to read data bytes.")]
    ReadDataError,
    #[error("Failed to read parity bytes.")]
    ReadParityError,
    #[error(transparent)]
    RSDecodeError(#[from] RSDecodeError),
    #[error(transparent)]
    TryFromIntError(#[from] TryFromIntError),
}
