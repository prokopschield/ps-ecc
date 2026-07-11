use crate::{Codeword, Cow};

impl<'lt> From<Cow<'lt>> for Codeword<'lt> {
    fn from(value: Cow<'lt>) -> Self {
        let range = 0..value.len();
        let codeword = value;

        Self { codeword, range }
    }
}
