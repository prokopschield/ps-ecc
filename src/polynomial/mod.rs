mod constants;
mod implementations;
mod methods;

use constants::POLYNOMIAL_MAX_COEFFICIENTS;

/// A polynomial over GF(256) with at most 255 coefficients.
///
/// Coefficients are stored from degree 0 upward in a fixed stack-allocated
/// array, so the type is `Copy` and never allocates.
#[derive(Clone, Copy, Eq)]
pub struct Polynomial {
    coefficients: [u8; POLYNOMIAL_MAX_COEFFICIENTS],
    degree: u8,
}
