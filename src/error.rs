#![allow(clippy::module_name_repetitions)]

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
