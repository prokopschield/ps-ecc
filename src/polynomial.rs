use ps_buffer::{Buffer, ToBuffer};

use crate::{
    error::PolynomialError,
    finite_field::{add, div, mul, sub},
};

/// Multiplies two polynomials.
pub fn poly_mul(p1: &[u8], p2: &[u8]) -> Result<Buffer, PolynomialError> {
    let mut result = Buffer::alloc(p1.len() + p2.len() - 1)?;
    for i in 0..p1.len() {
        for j in 0..p2.len() {
            result[i + j] = add(result[i + j], mul(p1[i], p2[j]));
        }
    }
    trim_leading_zeros(&mut result);
    Ok(result)
}

/// Evaluates a polynomial at a given point.
pub fn poly_eval(poly: &[u8], x: u8) -> u8 {
    let mut result = 0u8;
    for &coef in poly.iter().rev() {
        result = add(mul(result, x), coef);
    }
    result
}

/// Evaluates a polynomial split across two buffers at a given point.
/// The polynomial is conceptually [`low_degree_terms` || `high_degree_terms`].
/// This is equivalent to [`poly_eval`] but for when coefficients are stored
/// in separate buffers (e.g., parity bytes and message bytes).
pub fn poly_eval_detached(low_degree_terms: &[u8], high_degree_terms: &[u8], x: u8) -> u8 {
    let mut result = 0u8;

    for &coef in high_degree_terms.iter().rev() {
        result = add(mul(result, x), coef);
    }

    for &coef in low_degree_terms.iter().rev() {
        result = add(mul(result, x), coef);
    }

    result
}

/// Computes the remainder of polynomial division.
pub fn poly_rem(dividend: Buffer, divisor: &[u8]) -> Result<Buffer, PolynomialError> {
    use PolynomialError::ZeroDivisor;
    let mut rem = dividend;
    let divisor_deg = degree(divisor).ok_or(ZeroDivisor)?;
    while let Some(deg) = degree(&rem) {
        if deg < divisor_deg {
            break;
        }
        let lead_coef = rem[deg];
        let shift = deg - divisor_deg;
        for (i, item) in divisor.iter().enumerate() {
            let idx = shift + i;
            if idx < rem.len() {
                rem[idx] = sub(rem[idx], mul(lead_coef, *item));
            }
        }
        trim_leading_zeros(&mut rem);
    }
    Ok(rem)
}

/// Polynomial division returning quotient and remainder.
pub fn poly_div(dividend: &[u8], divisor: &[u8]) -> Result<(Buffer, Buffer), PolynomialError> {
    use PolynomialError::ZeroDivisor;

    let divisor_deg = degree(divisor).ok_or(ZeroDivisor)?;
    let divisor_lc = divisor[divisor_deg];

    let dividend_deg = degree(dividend).unwrap_or(0);
    let max_quot_deg = if dividend_deg >= divisor_deg {
        dividend_deg - divisor_deg
    } else {
        0
    };
    let mut quot = Buffer::alloc(max_quot_deg + 1)?;
    let mut rem = dividend.to_buffer()?;

    while let Some(deg) = degree(&rem) {
        if deg < divisor_deg {
            break;
        }
        let lead_coef = rem[deg];
        let ratio = div(lead_coef, divisor_lc)?;
        let shift = deg - divisor_deg;

        // Assert that shift is within bounds
        debug_assert!(shift < quot.len(), "Quotient index out of bounds");

        quot[shift] = ratio;

        // Subtract the scaled divisor from remainder
        for (i, &divisor_coef) in divisor.iter().enumerate() {
            let idx = shift + i;
            if idx < rem.len() {
                rem[idx] = sub(rem[idx], mul(ratio, divisor_coef));
            }
        }
        trim_leading_zeros(&mut rem);
    }

    trim_leading_zeros(&mut quot);
    Ok((quot, rem))
}

/// Evaluates the derivative of a polynomial at a point (in characteristic 2).
pub fn poly_eval_deriv(poly: &[u8], x: u8) -> u8 {
    let mut result = 0u8;
    let x_sq = mul(x, x);
    let mut x_pow = 1u8;
    for k in 0..=((poly.len() - 1) / 2) {
        let idx = 2 * k + 1;
        if idx < poly.len() {
            result = add(result, mul(poly[idx], x_pow));
            x_pow = mul(x_pow, x_sq);
        }
    }
    result
}

/// Subtracts two polynomials (same as addition in GF(2)).
pub fn poly_sub(p1: &[u8], p2: &[u8]) -> Result<Buffer, PolynomialError> {
    let len = p1.len().max(p2.len());
    let mut result = Buffer::alloc(len)?;

    for (i, coef) in result.iter_mut().enumerate() {
        let a = p1.get(i).copied().unwrap_or(0);
        let b = p2.get(i).copied().unwrap_or(0);
        *coef = add(a, b);
    }

    trim_leading_zeros(&mut result);
    Ok(result)
}

/// Computes the degree of a polynomial.
fn degree(poly: &[u8]) -> Option<usize> {
    poly.iter().rposition(|&x| x != 0)
}

/// Removes leading zeros from a polynomial.
fn trim_leading_zeros(poly: &mut Buffer) {
    poly.truncate(degree(poly).unwrap_or(0).saturating_add(1));
}

#[cfg(test)]
mod tests {
    use crate::PolynomialError;

    use super::poly_div;

    #[test]
    fn try_poly_div() -> Result<(), PolynomialError> {
        let (q, r) = poly_div(&[1, 2, 3, 4, 5], &[1, 1, 0, 0])?;

        assert_eq!(q.slice(..), &[0, 2, 1, 5]);
        assert_eq!(r.slice(..), &[1]);

        Ok(())
    }
}
