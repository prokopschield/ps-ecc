//! Optimized Extended Euclidean Algorithm for Reed-Solomon decoding.
//!
//! This implementation minimizes memory allocations by using stack-allocated
//! polynomials and ring buffer index swapping instead of data movement.

use crate::{EuclideanError, Polynomial};

/// Extended Euclidean algorithm for Reed-Solomon decoding.
///
/// Given syndrome coefficients and parity count `t`, computes the error
/// locator polynomial (sigma) and error evaluator polynomial (omega).
///
/// # Arguments
///
/// * `syndromes` - Syndrome coefficients in ascending degree order
/// * `t` - Number of parity symbols (must satisfy `2*t <= 255`)
///
/// # Returns
///
/// A tuple `(sigma, omega)` where:
/// * `sigma` is the error locator polynomial
/// * `omega` is the error evaluator polynomial
///
/// # Errors
///
/// Returns an error if:
/// * `syndromes` exceeds 255 bytes
/// * Division by zero occurs (indicates invalid input)
pub fn euclidean(syndromes: &[u8], t: u8) -> Result<(Polynomial, Polynomial), EuclideanError> {
    if syndromes.len() > 255 {
        return Err(EuclideanError::SyndromesTooLong(syndromes.len()));
    }

    // Ring buffers: r[idx] is current, r[idx ^ 1] is previous
    let mut r: [Polynomial; 2] = [Polynomial::default(), Polynomial::default()];
    let mut t_poly: [Polynomial; 2] = [Polynomial::default(), Polynomial::default()];

    // r[0] = x^(2t)
    let two_t = usize::from(t) * 2;

    #[allow(clippy::cast_possible_truncation)]
    r[0].set(two_t as u8, 1);

    // r[1] = syndrome polynomial
    r[1] = Polynomial::try_from(syndromes)?;

    // t[0] = 0 (already default)
    // t[1] = 1
    t_poly[1].set(0, 1);

    // Reusable quotient polynomial
    let mut q = Polynomial::default();

    // idx points to current r and t
    let mut idx: usize = 1;

    while r[idx].degree() >= t {
        // Split arrays to satisfy borrow checker
        // idx is always 0 or 1, so one of (r0, r1) is (prev, curr)
        let (r_left, r_right) = r.split_at_mut(1);
        let (r_prev, r_curr) = if idx == 1 {
            (&mut r_left[0], &r_right[0])
        } else {
            (&mut r_right[0], &r_left[0])
        };

        let (t_left, t_right) = t_poly.split_at_mut(1);
        let (t_prev, t_curr) = if idx == 1 {
            (&mut t_left[0], &t_right[0])
        } else {
            (&mut t_right[0], &t_left[0])
        };

        // r[prev] <- r[prev] mod r[curr], q <- quotient
        r_prev.div_rem_inplace(r_curr, &mut q)?;

        // t[prev] <- t[prev] ^ (q * t[curr])
        t_prev.mul_xor_assign(&q, t_curr);

        // Flip index: what was prev becomes curr
        idx ^= 1;
    }

    // Return (sigma, omega) = (t[idx], r[idx])
    Ok((t_poly[idx], r[idx]))
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use super::*;
    use crate::finite_field::{inv, mul, ANTILOG_TABLE};

    #[test]
    fn euclidean_zero_syndromes() {
        // All-zero syndromes indicate no errors
        let syndromes = [0u8; 4];
        let t = 2;

        let (sigma, omega) = euclidean(&syndromes, t).expect("should succeed");

        // With no errors, sigma should be constant (degree 0)
        assert_eq!(sigma.degree(), 0);
        // omega should be zero polynomial
        assert_eq!(omega.coefficients(), &[0]);
    }

    #[test]
    fn euclidean_single_error_sigma_has_one_root() {
        // For a single error at position j, sigma(x) = 1 - alpha^j * x
        // The syndrome for a single error e at position j is: S_i = e * alpha^(i*j)
        // Let's use error value e=1 at position j=0, so S_i = alpha^0 = 1 for all i
        let syndromes = [1u8; 4]; // Single error at position 0 with value 1
        let t = 2;

        let (sigma, _omega) = euclidean(&syndromes, t).expect("should succeed");

        // sigma should have degree 1 (one error)
        assert!(sigma.degree() <= 1);

        // Verify sigma(alpha^0) = 0 (error at position 0)
        // alpha^(255-0) mod 255 = alpha^0 = 1
        let x = ANTILOG_TABLE[255 % 255].get(); // alpha^0 = 1

        let sigma_at_x = Polynomial::eval_coefficients_at(sigma.coefficients(), x);

        // After normalization, sigma(x) should be zero at the error position
        // Note: euclidean returns unnormalized sigma, so we normalize first
        let scale = inv(sigma[0u8]).expect("sigma[0] should be non-zero").get();
        let _normalized_eval = mul(sigma_at_x, scale);

        // Due to how the algorithm works, the root check is on sigma directly
        // sigma(alpha^(-j)) = 0 for error at position j
        assert_eq!(
            Polynomial::eval_coefficients_at(sigma.coefficients(), 1),
            0,
            "sigma should have a root at alpha^0 for error at position 0"
        );
    }

    #[test]
    fn euclidean_verifies_degree_bounds() {
        // The extended Euclidean algorithm guarantees:
        // - deg(sigma) <= t (error locator has at most t roots)
        // - deg(omega) < t (loop exits when remainder degree < t)
        let syndromes = [3u8, 5, 7, 11]; // Arbitrary non-zero syndromes
        let t = 2;

        let (sigma, omega) = euclidean(&syndromes, t).expect("should succeed");

        // Verify deg(sigma) <= t (can correct at most t errors)
        assert!(
            sigma.degree() <= t,
            "sigma degree {} exceeds t={}",
            sigma.degree(),
            t
        );

        // Verify deg(omega) < t (algorithm invariant)
        assert!(
            omega.degree() < t || omega.coefficients() == [0],
            "omega degree {} should be < t={}",
            omega.degree(),
            t
        );
    }

    #[test]
    fn euclidean_with_t_equals_one() {
        // Minimum parity: can correct 1 error
        let syndromes = [42u8, 0]; // 2t = 2 syndromes
        let t = 1;

        let (sigma, omega) = euclidean(&syndromes, t).expect("should succeed");

        assert!(sigma.degree() <= t);
        assert!(omega.degree() < t || omega.coefficients() == [0]);
    }

    #[test]
    fn euclidean_with_large_t() {
        // Larger parity count
        let syndromes = [1u8; 20]; // t = 10
        let t = 10;

        let (sigma, omega) = euclidean(&syndromes, t).expect("should succeed");

        assert!(sigma.degree() <= t);
        assert!(omega.degree() < t || omega.coefficients() == [0]);
    }

    #[test]
    fn euclidean_minimal_valid_case() {
        // Minimal meaningful case: t=1, 2 syndromes
        let syndromes = [42u8, 17];
        let t = 1;

        let (sigma, omega) = euclidean(&syndromes, t).expect("should succeed");

        assert!(sigma.degree() <= t);
        assert!(omega.degree() < t || omega.coefficients() == [0]);
    }

    #[test]
    fn euclidean_syndrome_with_trailing_zeros() {
        // Syndromes with trailing zeros (leading zeros in polynomial sense)
        let syndromes = [1u8, 2, 0, 0]; // Effective degree 1, but length 4
        let t = 2;

        let (sigma, omega) = euclidean(&syndromes, t).expect("should succeed");

        // sigma degree bounded by t
        assert!(sigma.degree() <= t);

        // omega degree bounded by t (loop exits when deg(r) < t)
        assert!(
            omega.degree() < t || omega.coefficients() == [0],
            "omega degree {} should be < t={}",
            omega.degree(),
            t
        );
    }

    #[test]
    fn euclidean_too_long() {
        let syndromes = [0u8; 256];

        let result = euclidean(&syndromes, 1);

        assert!(matches!(result, Err(EuclideanError::SyndromesTooLong(256))));
    }

    #[test]
    fn euclidean_max_valid_length() {
        // 255 bytes is the maximum valid length
        let syndromes = [0u8; 255];
        let t = 127;

        let result = euclidean(&syndromes, t);

        assert!(result.is_ok());
    }

    #[test]
    fn euclidean_output_is_deterministic() {
        // Same input should always produce same output
        let syndromes = [7u8, 13, 19, 23, 29, 31];
        let t = 3;

        let (sigma1, omega1) = euclidean(&syndromes, t).expect("should succeed");
        let (sigma2, omega2) = euclidean(&syndromes, t).expect("should succeed");

        assert_eq!(sigma1, sigma2);
        assert_eq!(omega1, omega2);
    }

    #[test]
    fn euclidean_integration_with_rs_syndromes() {
        // Use actual RS syndrome computation for a realistic test
        // Create a simple codeword with a known error pattern
        use crate::ReedSolomon;

        let rs = ReedSolomon::new(2).expect("valid RS");

        // Encode a simple message
        let message = b"test";
        let mut encoded = rs.encode(message).expect("encoding succeeds");

        // Introduce a single error
        encoded[0] ^= 0x42; // Flip some bits in first byte

        // Compute syndromes
        let syndromes =
            ReedSolomon::compute_syndromes(rs.parity_bytes(), &encoded).expect("syndromes");

        // Run euclidean algorithm
        let t = rs.parity_bytes() / 2;
        let (sigma, omega) = euclidean(&syndromes, t).expect("euclidean succeeds");

        // Verify sigma has at most t roots
        assert!(sigma.degree() <= t, "sigma degree should be <= t");

        // Verify omega degree constraint (deg(omega) < t)
        assert!(
            omega.degree() < t || omega.coefficients() == [0],
            "omega degree {} should be < t={}",
            omega.degree(),
            t
        );
    }
}
