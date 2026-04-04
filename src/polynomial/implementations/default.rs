use crate::Polynomial;

impl Default for Polynomial {
    fn default() -> Self {
        Self {
            coefficients: [0; _],
            degree: 0,
        }
    }
}
