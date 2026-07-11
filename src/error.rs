#![allow(clippy::module_name_repetitions)]

use std::num::TryFromIntError;

use ps_buffer::BufferError;
use thiserror::Error;

/// Errors of GF(256) field arithmetic.
#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum GFError {
    /// The divisor of a field division was zero.
    #[error("Division by zero is undefined.")]
    DivByZero,
}

/// Errors returned when constructing a [`Polynomial`](crate::Polynomial)
/// from a slice.
#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialFromSliceError {
    /// The slice holds more coefficients than a polynomial can store.
    #[error("Slice was {size} bytes, max 255 allowed in GF(256).")]
    TooLong { size: usize },
}

/// Errors returned by
/// [`Polynomial::set_coefficients`](crate::Polynomial::set_coefficients).
#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialSetCoefficientsError {
    /// The coefficient range extends past the maximum coefficient index.
    #[error("Range {offset}..{end} exceeds maximum coefficient index 254.")]
    OutOfBounds { offset: u8, end: usize },
}

/// Errors of polynomial multiplication.
#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialMulError {
    /// The degree of the product exceeds
    /// [`Polynomial::MAX_DEGREE`](crate::Polynomial::MAX_DEGREE).
    #[error("Result degree exceeds maximum of 254.")]
    DegreeOverflow,
}

/// Errors of polynomial division and remainder.
#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialDivError {
    /// The divisor was the zero polynomial.
    #[error("Division by zero polynomial.")]
    ZeroDivisor,
    /// Propagated from GF(256) arithmetic.
    #[error(transparent)]
    GFError(#[from] GFError),
}

/// Errors of polynomial XOR with a coefficient iterator.
#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolynomialXorError {
    /// The iterator yielded more coefficients than a polynomial can store.
    #[error("Coefficient iterator exceeded maximum length of 255.")]
    TooManyCoefficients,
}

/// Errors returned by [`euclidean`](crate::euclidean).
#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum EuclideanError {
    /// The error-correction capability `t` requires the polynomial `x^(2t)`,
    /// whose degree exceeds [`Polynomial::MAX_DEGREE`](crate::Polynomial::MAX_DEGREE).
    #[error("Error-correction capability {t} exceeds the maximum of 127.")]
    CapabilityTooHigh { t: u8 },
    /// Propagated from polynomial division.
    #[error(transparent)]
    PolynomialDiv(#[from] PolynomialDivError),
    /// Propagated from polynomial multiplication.
    #[error(transparent)]
    PolynomialMul(#[from] PolynomialMulError),
}

/// Errors returned by [`ReedSolomon::new`](crate::ReedSolomon::new) and by
/// the detached-parity methods when the parity slice itself is invalid.
#[derive(Error, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum RSConstructorError {
    /// A detached parity slice holds an odd number of bytes; parity always
    /// comprises two bytes per correctable error.
    #[error("Parity byte count {0} is odd; parity comprises two bytes per correctable error.")]
    OddParityLength(usize),
    /// The requested error-correction capability exceeds
    /// [`MAX_PARITY`](crate::MAX_PARITY).
    #[error("Parity count must be <= 63.")]
    ParityTooHigh,
}

/// Errors returned by
/// [`ReedSolomon::generate_parity`](crate::ReedSolomon::generate_parity).
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSGenerateParityError {
    /// Propagated from polynomial division.
    #[error(transparent)]
    Division(#[from] PolynomialDivError),
    /// Propagated from assembling the message polynomial; the message and
    /// parity together exceed a single codeword.
    #[error(transparent)]
    SetCoefficients(#[from] PolynomialSetCoefficientsError),
}

/// Errors returned by [`ReedSolomon::encode`](crate::ReedSolomon::encode).
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSEncodeError {
    /// Propagated from buffer allocation.
    #[error(transparent)]
    BufferError(#[from] BufferError),
    /// Propagated from parity generation.
    #[error(transparent)]
    RSGenerateParityError(#[from] RSGenerateParityError),
}

/// Errors returned by
/// [`ReedSolomon::compute_errors`](crate::ReedSolomon::compute_errors).
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSComputeErrorsError {
    /// Propagated from GF(256) arithmetic.
    #[error(transparent)]
    GFError(#[from] GFError),
    /// Propagated from the extended Euclidean algorithm.
    #[error(transparent)]
    EuclideanError(#[from] EuclideanError),
    /// The codeword holds more errors than the parity can correct.
    #[error("Too many errors; input is unrecoverable.")]
    TooManyErrors,
    /// Forney's formula would divide by zero; not expected to occur.
    #[error("The error locator derivative evaluated to zero.")]
    ZeroErrorLocatorDerivative,
}

/// Errors returned by the decoding and correction methods of
/// [`ReedSolomon`](crate::ReedSolomon).
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum RSDecodeError {
    /// Propagated from buffer allocation.
    #[error(transparent)]
    BufferError(#[from] BufferError),
    /// The input holds fewer bytes than the parity, so it cannot be a
    /// codeword.
    #[error("Input length {received} is less than the parity length {parity_bytes}.")]
    InsufficientLength { parity_bytes: u8, received: usize },
    /// Propagated from error computation.
    #[error(transparent)]
    RSComputeErrorsError(#[from] RSComputeErrorsError),
    /// Propagated from parity-slice validation.
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    /// The corrected bytes still fail validation.
    #[error("Too many errors to correct. Error computation nevertheless returned a valid polynomial, which is unlikely. Usually you'll get RSComputeErrorsError(TooManyErrors) instead.")]
    TooManyErrors,
    /// The input holds more bytes than a single codeword can carry.
    #[error(transparent)]
    TryFromIntError(#[from] TryFromIntError),
}

/// Errors returned by the free [`encode`](crate::encode) function.
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum EncodeError {
    /// Propagated from long ECC encoding.
    #[error(transparent)]
    LongEccEncodeError(#[from] LongEccEncodeError),
    /// Propagated from codec construction.
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    /// Propagated from single-codeword encoding.
    #[error(transparent)]
    RSEncodeError(#[from] RSEncodeError),
}

/// Errors returned by the free [`decode`](crate::decode) function.
#[derive(Error, Debug, Clone)]
pub enum DecodeError {
    /// The input is too short to carry the requested parity.
    #[error("Insufficient input bytes for parity count of {0}: {0} * 2 > {1}.")]
    InsufficientParityBytes(u8, u8),
    /// Propagated from long ECC decoding.
    #[error(transparent)]
    LongEccDecodeError(#[from] LongEccDecodeError),
    /// Propagated from codec construction.
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    /// Propagated from single-codeword decoding.
    #[error(transparent)]
    RSDecodeError(#[from] RSDecodeError),
}

/// Umbrella error over encoding and decoding.
#[derive(Error, Debug, Clone)]
pub enum EccError {
    /// Propagated from encoding.
    #[error(transparent)]
    EncodeError(#[from] EncodeError),
    /// Propagated from decoding.
    #[error(transparent)]
    DecodeError(#[from] DecodeError),
}

/// Errors returned by [`LongEccHeader::new`](crate::LongEccHeader::new).
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccHeaderConstructorError {
    /// Propagated from generating the header parity.
    #[error("Generating parity failed: {0}")]
    GenerateParity(#[from] RSGenerateParityError),
    /// The full length does not match the codeword length derived from the
    /// message length and segment geometry.
    #[error("Full length {0} does not match the derived codeword length {1}.")]
    InvalidFullLength(u32, u64),
    /// The header and message do not fit within the full length.
    #[error("Message length {0} does not fit within full length {1}.")]
    InvalidMessageLength(u32, u32),
    /// The error-correction capability exceeds
    /// [`MAX_PARITY`](crate::MAX_PARITY).
    #[error("Invalid parity count: {0}.")]
    InvalidParityCount(u8),
    /// The parity bytes leave no room for new data within a segment.
    #[error("Invalid segment-to-parity ratio: {0} <= 2 * {1}.")]
    InvalidSegmentParityRatio(u8, u8),
}

/// Errors returned by
/// [`LongEccHeader::from_bytes`](crate::LongEccHeader::from_bytes).
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccHeaderFromBytesError {
    /// The magic number does not identify a long ECC header.
    #[error("Incorrect magic number: {0:x}.")]
    IncorrectMagic(u16),

    /// The encoding version is not supported.
    #[error("Incorrect version number: {0}.")]
    InvalidVersion(u8),

    /// The header checksum mismatches even after error correction.
    #[error("Header checksum incorrect.")]
    IncorrectChecksum,

    /// Propagated from correcting the header bytes with the header parity.
    #[error("Header error correction failed: {0}")]
    CorrectionFailed(#[from] RSDecodeError),

    /// The header and message do not fit within the full length.
    #[error("Message length {0} does not fit within full length {1}.")]
    InvalidMessageLength(u32, u32),

    /// The parity bytes leave no room for new data within a segment.
    #[error("Invalid segment-to-parity ratio: {0} <= 2 * {1}.")]
    InvalidSegmentParityRatio(u8, u8),

    /// The full length does not match the codeword length derived from the
    /// message length and segment geometry.
    #[error("Full length {0} does not match the derived codeword length {1}.")]
    InvalidFullLength(u32, u64),
}

/// Errors returned by
/// [`LongEccHeader::from_byte_slice`](crate::LongEccHeader::from_byte_slice).
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccHeaderFromByteSliceError {
    /// The slice holds fewer than the 32 bytes a header occupies.
    #[error("Insufficient bytes for header: got {0}, need 32.")]
    InsufficientBytes(usize),

    /// Propagated from parsing the 32-byte header.
    #[error(transparent)]
    FromBytes(#[from] LongEccHeaderFromBytesError),
}

/// Errors of long ECC encoding.
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum LongEccEncodeError {
    /// Propagated from buffer allocation.
    #[error(transparent)]
    BufferError(#[from] BufferError),
    /// The requested error-correction capability exceeds
    /// [`MAX_PARITY`](crate::MAX_PARITY).
    #[error("Parity {0} >= 64, which is too high.")]
    InvalidParity(u8),
    /// The parity bytes leave no room for new data within a segment.
    #[error("Invalid segment-to-parity ratio: {0} <= 2 * {1}.")]
    InvalidSegmentParityRatio(u8, u8),
    /// Propagated from header construction.
    #[error("Long ECC header construction failed: {0}")]
    LongEccHeaderCtor(#[from] LongEccHeaderConstructorError),
    /// Propagated from codec construction.
    #[error(transparent)]
    RSConstructorError(#[from] RSConstructorError),
    /// Propagated from parity generation.
    #[error(transparent)]
    RSGenerateParityError(#[from] RSGenerateParityError),
    /// The encoded codeword would exceed `u32::MAX` bytes.
    #[error(transparent)]
    TryFromIntError(#[from] TryFromIntError),
}

/// Errors of long ECC decoding.
#[derive(Error, Debug, Clone)]
pub enum LongEccDecodeError {
    /// Propagated from buffer allocation.
    #[error(transparent)]
    BufferError(#[from] BufferError),
    /// The checksum still mismatches after error correction.
    #[error("Integrity check failed after correction.")]
    IntegrityCheckFailed,
    /// The codeword is structurally invalid.
    #[error("Codeword is invalid.")]
    InvalidCodeword,
    /// Propagated from header parsing.
    #[error("Failed to decode header: {0}")]
    HeaderDecode(#[from] LongEccHeaderFromByteSliceError),
    /// The data bytes recorded in the header lie outside the buffer.
    #[error("Failed to read data bytes.")]
    ReadDataError,
    /// The parity bytes recorded in the header lie outside the buffer.
    #[error("Failed to read parity bytes.")]
    ReadParityError,
    /// Propagated from single-codeword decoding of a segment.
    #[error(transparent)]
    RSDecodeError(#[from] RSDecodeError),
    /// A header length field does not fit the platform's `usize`.
    #[error(transparent)]
    TryFromIntError(#[from] TryFromIntError),
}
