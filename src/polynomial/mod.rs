mod constants;
mod implementations;
mod methods;

use constants::POLYNOMIAL_MAX_COEFFICIENTS;
pub use methods::{poly_mul, poly_rem};

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct Polynomial {
    coefficients: [u8; POLYNOMIAL_MAX_COEFFICIENTS],
    degree: u8,
}
