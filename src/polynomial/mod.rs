mod constants;
mod implementations;
mod methods;

use constants::POLYNOMIAL_MAX_COEFFICIENTS;
pub use methods::{poly_div, poly_eval_deriv, poly_eval_detached, poly_mul, poly_rem, poly_sub};

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct Polynomial {
    coefficients: [u8; POLYNOMIAL_MAX_COEFFICIENTS],
    degree: u8,
}
