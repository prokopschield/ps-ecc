//! Reed-Solomon error-correcting codes over GF(256).
//!
//! The crate offers three layers of API:
//!
//! - The free functions [`encode`], [`decode`], and [`validate`] pick the
//!   appropriate format automatically: messages that fit a single
//!   Reed-Solomon codeword (at most `255 - 2 * parity` bytes) are encoded
//!   directly, and longer messages use the long ECC format, which splits
//!   the message into overlapping 255-byte segments behind a checksummed
//!   header.
//! - [`ReedSolomon`] encodes, validates, and corrects single codewords,
//!   either attached (`parity || data` in one slice) or detached (parity
//!   and data as separate slices).
//! - [`Polynomial`] and [`euclidean`] expose the underlying GF(256)
//!   polynomial arithmetic.
//!
//! The `parity` parameter throughout is the error-correction capability:
//! the number of byte errors a codeword can recover from, each costing two
//! parity bytes. It is capped at [`MAX_PARITY`] (63).

mod codeword;
mod cow;
mod error;
mod euclidean;
mod finite_field;
mod long;
mod methods;
mod polynomial;
mod reed_solomon;

pub use codeword::Codeword;
pub use cow::Cow;
pub use error::*;
pub use euclidean::euclidean;
pub use long::{
    LongEccHeader, LongEccHeaderConstructorError, LongEccHeaderFromByteSliceError,
    LongEccHeaderFromBytesError, OverlapFactor,
};
pub use methods::*;
pub use polynomial::Polynomial;
pub use reed_solomon::*;
