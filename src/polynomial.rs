use ps_buffer::{Buffer, ToBuffer};

use crate::{
    error::PolynomialError,
    finite_field::{add, div, mul},
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

pub fn poly_eval_detached(left: &[u8], right: &[u8], x: u8) -> u8 {
    let mut result = 0u8;
    for &coef in right.iter().rev() {
        result = add(mul(result, x), coef);
    }
    for &coef in left.iter().rev() {
        result = add(mul(result, x), coef);
    }
    result
}

/// Computes the remainder of polynomial division.
pub fn poly_rem(dividend: &[u8], divisor: &[u8]) -> Result<Buffer, PolynomialError> {
    use PolynomialError::ZeroDivisor;
    let mut rem = dividend.to_buffer()?;
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
                rem[idx] = add(rem[idx], mul(lead_coef, *item));
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
    let mut quot = Buffer::alloc(dividend.len().saturating_sub(divisor.len()) + 1)?;
    let mut rem = dividend.to_buffer()?;
    while let Some(deg) = degree(&rem) {
        if deg < divisor_deg {
            break;
        }
        let lead_coef = rem[deg];
        let ratio = div(lead_coef, divisor_lc)?;
        let shift = deg - divisor_deg;
        quot[shift] = ratio;
        for (i, item) in divisor.iter().enumerate() {
            let idx = shift + i;
            if idx < rem.len() {
                rem[idx] = add(rem[idx], mul(ratio, *item));
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
#[allow(clippy::needless_range_loop)]
pub fn poly_sub(p1: &[u8], p2: &[u8]) -> Result<Buffer, PolynomialError> {
    let len = p1.len().max(p2.len());
    let mut result = Buffer::alloc(len)?;
    for i in 0..len {
        let a = p1.get(i).copied().unwrap_or(0);
        let b = p2.get(i).copied().unwrap_or(0);
        result[i] = add(a, b);
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
    poly.truncate(degree(poly).unwrap_or(1).saturating_add(1));
}
