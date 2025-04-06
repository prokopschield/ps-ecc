#![allow(clippy::module_name_repetitions)]

use std::num::TryFromIntError;

use ps_buffer::BufferError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum GFError {
    #[error("Division by zero is undefined.")]
    DivByZero,
}

#[derive(Error, Debug)]
pub enum PolynomialError {
    #[error(transparent)]
    GFError(#[from] GFError),
    #[error("Divisor cannot be zero.")]
    ZeroDivisor,
}

#[derive(Error, Debug)]
pub enum RSConstructorError {
    #[error("Parity count must be <= 127.")]
    ParityTooHigh,
}

#[derive(Error, Debug)]
pub enum RSEncodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    PolynomialError(#[from] PolynomialError),
}

#[derive(Error, Debug)]
pub enum RSDecodeError {
    #[error(transparent)]
    BufferError(#[from] BufferError),
    #[error(transparent)]
    GFError(#[from] GFError),
    #[error(transparent)]
    PolynomialError(#[from] PolynomialError),
    #[error("Too many errors to correct.")]
    TooManyErrors,
    #[error("Derivative evaluated to zero.")]
    ZeroDerivative,
}

#[derive(Error, Debug)]
pub enum EncodeError {
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    #[error(transparent)]
    RSEncodeError(#[from] RSEncodeError),
}

#[derive(Error, Debug)]
pub enum DecodeError {
    #[error("Input too large, {0} > 255.")]
    InputTooLarge(u32),
    #[error("Insufficient input bytes for parity count of {0}: {0} * 2 > {1}.")]
    InsufficientParityBytes(u8, u8),
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    #[error(transparent)]
    RSDecodeError(#[from] RSDecodeError),
    #[error(transparent)]
    TryFromIntError(#[from] TryFromIntError),
}
#[derive(Error, Debug)]
pub enum EccError {
    #[error(transparent)]
    EncodeError(#[from] EncodeError),
    #[error(transparent)]
    DecodeError(#[from] DecodeError),
}
