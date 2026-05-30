//! Compile-time generated Reed-Solomon generator polynomials.
//!
//! The generator polynomial for `t` correctable errors (2t parity bytes) is:
//! `g(x) = (x - α^1)(x - α^2)...(x - α^{2t})`
//!
//! Polynomials are stored in a flat buffer to eliminate padding waste.

use crate::{
    finite_field::{mul, ANTILOG_TABLE},
    MAX_PARITY,
};

/// Maximum parity as usize for array indexing.
const MAX_PARITY_USIZE: usize = MAX_PARITY as usize;

/// Total coefficients: Σ(2p + 1) for p = `0..=MAX_PARITY` = (`MAX_PARITY` + 1)² = 4096
const TOTAL_COEFFS: usize = (MAX_PARITY_USIZE + 1) * (MAX_PARITY_USIZE + 1);

/// Flat buffer containing all generator polynomials concatenated.
/// Polynomial for parity `p` starts at offset `p²` and has length `2p + 1`.
const GENERATOR_DATA: [u8; TOTAL_COEFFS] = compute_generator_data();

const fn compute_generator_data() -> [u8; TOTAL_COEFFS] {
    let mut data = [0u8; TOTAL_COEFFS];

    // Temporary workspace for building polynomials incrementally
    let mut poly = [0u8; 2 * MAX_PARITY_USIZE + 1];

    poly[0] = 1; // g_0(x) = 1

    // Copy g_0 to data
    data[0] = 1;

    let mut parity = 1;

    while parity <= MAX_PARITY_USIZE {
        // Multiply by (x + α^{2*parity-1}) and (x + α^{2*parity})
        let root1 = ANTILOG_TABLE[2 * parity - 1].get();
        let root2 = ANTILOG_TABLE[2 * parity].get();

        // Multiply by (x + root1)
        let deg1 = 2 * parity - 1;
        let mut j = deg1;

        while j > 0 {
            poly[j] = poly[j - 1] ^ mul(poly[j], root1);
            j -= 1;
        }

        poly[0] = mul(poly[0], root1);

        // Multiply by (x + root2)
        let deg2 = 2 * parity;

        j = deg2;

        while j > 0 {
            poly[j] = poly[j - 1] ^ mul(poly[j], root2);
            j -= 1;
        }

        poly[0] = mul(poly[0], root2);

        // Copy to flat buffer at correct offset
        let offset = parity * parity;
        let len = 2 * parity + 1;
        let mut i = 0;

        while i < len {
            data[offset + i] = poly[i];
            i += 1;
        }

        parity += 1;
    }

    data
}

/// Returns the generator polynomial for the given parity (error correction capability).
///
/// # Panics
///
/// Panics if `parity > 63`.
#[inline]
pub(super) fn generator_poly(parity: u8) -> &'static [u8] {
    let p = parity as usize;
    let start = p * p;
    let len = 2 * p + 1;

    &GENERATOR_DATA[start..start + len]
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation)]
mod tests {
    use super::*;

    #[test]
    fn generator_poly_parity_zero() {
        assert_eq!(generator_poly(0), &[1]);
    }

    #[test]
    fn generator_poly_parity_one() {
        // g(x) = (x + α^1)(x + α^2) = x² + (α^1 ⊕ α^2)x + α^3
        let g = generator_poly(1);

        assert_eq!(g.len(), 3);

        let alpha1 = ANTILOG_TABLE[1].get();
        let alpha2 = ANTILOG_TABLE[2].get();
        let alpha3 = ANTILOG_TABLE[3].get();

        assert_eq!(g[0], alpha3);
        assert_eq!(g[1], alpha1 ^ alpha2);
        assert_eq!(g[2], 1);
    }

    #[test]
    fn generator_poly_is_monic() {
        for p in 0..=MAX_PARITY_USIZE {
            let g = generator_poly(p as u8);
            let degree = 2 * p;

            assert_eq!(g[degree], 1, "parity {p} should be monic");
        }
    }

    #[test]
    fn generator_poly_correct_length() {
        for p in 0..=MAX_PARITY_USIZE {
            let expected = 2 * p + 1;

            assert_eq!(generator_poly(p as u8).len(), expected);
        }
    }

    #[test]
    fn generator_poly_has_correct_roots() {
        for p in 1..=8usize {
            let g = generator_poly(p as u8);

            for (i, alpha_i) in ANTILOG_TABLE.iter().take(2 * p + 1).enumerate().skip(1) {
                let eval = eval_poly(g, alpha_i.get());

                assert_eq!(eval, 0, "parity {p}: g(α^{i}) should be 0");
            }
        }
    }

    fn eval_poly(coeffs: &[u8], x: u8) -> u8 {
        let mut result = 0u8;

        for &c in coeffs.iter().rev() {
            result = mul(result, x) ^ c;
        }

        result
    }
}
