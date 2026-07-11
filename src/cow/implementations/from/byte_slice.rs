use crate::Cow;

impl<'lt> From<&'lt [u8]> for Cow<'lt> {
    fn from(value: &'lt [u8]) -> Self {
        Self::Borrowed(value)
    }
}
