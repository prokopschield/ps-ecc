mod constants;
mod generator;
mod methods;
mod types;

pub use constants::*;
pub use types::*;

/// A Reed-Solomon codec over GF(256) for codewords of at most 255 bytes.
///
/// The codec is parameterized by its error-correction capability: the
/// number of byte errors it can correct, each costing two parity bytes.
/// Codewords are laid out as `parity || data`.
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReedSolomon {
    parity: u8,
}
