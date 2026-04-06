mod coefficients;
mod degree;
mod div_rem;
mod div_rem_inplace;
mod eval_at;
mod eval_coefficient_iter_at;
mod eval_coefficient_slices_at;
mod eval_coefficients_at;
mod eval_coefficients_derivative_at;
mod eval_derivative_at;
mod from_slice;
mod iter;
mod mul_xor_assign;
mod set;
mod set_coefficients;
mod trim_degree;

use ps_buffer::Buffer;

use crate::{
    error::PolynomialError,
    finite_field::{add, mul, sub},
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

/// Computes the degree of a polynomial.
fn degree(poly: &[u8]) -> Option<usize> {
    poly.iter().rposition(|&x| x != 0)
}

/// Removes leading zeros from a polynomial.
fn trim_leading_zeros(poly: &mut Buffer) {
    poly.truncate(degree(poly).unwrap_or(0).saturating_add(1));
}
