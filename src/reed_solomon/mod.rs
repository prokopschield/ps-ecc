mod constants;
mod generator;
mod methods;
mod types;

pub use constants::*;
pub use types::*;

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReedSolomon {
    parity: u8,
}
