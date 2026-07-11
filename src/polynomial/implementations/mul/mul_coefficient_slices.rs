use crate::{
    error::PolynomialMulError,
    finite_field::{add, mul},
    Polynomial,
};

/// Core multiplication on trimmed coefficient slices.
///
/// Both slices must represent trimmed polynomials (no trailing zeros beyond degree).
/// Empty slices are treated as the zero polynomial.
pub(super) fn mul_coefficient_slices(
    lhs: &[u8],
    rhs: &[u8],
) -> Result<Polynomial, PolynomialMulError> {
    if lhs.is_empty() || rhs.is_empty() {
        return Ok(Polynomial::default());
    }

    let lhs_degree = lhs.len() - 1;
    let rhs_degree = rhs.len() - 1;

    let result_degree = lhs_degree + rhs_degree;

    if result_degree > Polynomial::MAX_DEGREE as usize {
        return Err(PolynomialMulError::DegreeOverflow);
    }

    #[allow(clippy::cast_possible_truncation)]
    let mut result = Polynomial {
        degree: result_degree as u8,
        ..Polynomial::default()
    };

    for (i, &lhs_coef) in lhs.iter().enumerate() {
        for (j, &rhs_coef) in rhs.iter().enumerate() {
            let k = i + j;

            result.coefficients[k] = add(result.coefficients[k], mul(lhs_coef, rhs_coef));
        }
    }

    result.trim_degree();

    Ok(result)
}
