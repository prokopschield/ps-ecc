use crate::{long::header::magic::LONG_ECC_HEADER_PARITY, ReedSolomon};

pub const RS: ReedSolomon = match ReedSolomon::new(LONG_ECC_HEADER_PARITY) {
    Ok(codec) => codec,
    Err(_) => panic!("Failed to initialize ReedSolomon codec at compile time."),
};
