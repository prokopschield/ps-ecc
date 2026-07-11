use crate::Cow;

impl PartialEq for Cow<'_> {
    fn eq(&self, other: &Self) -> bool {
        **self == **other
    }
}

impl Eq for Cow<'_> {}

#[cfg(test)]
mod tests {
    use ps_buffer::Buffer;

    use crate::Cow;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_eq_across_variants_with_same_content() -> Result<(), TestError> {
        let borrowed = Cow::Borrowed(b"abc");
        let owned = Cow::Owned(Buffer::from_slice(b"abc")?.share());

        assert_eq!(borrowed, owned);

        Ok(())
    }

    #[test]
    fn test_ne_for_different_content() -> Result<(), TestError> {
        let borrowed = Cow::Borrowed(b"abc");
        let owned = Cow::Owned(Buffer::from_slice(b"abd")?.share());

        assert_ne!(borrowed, owned);

        Ok(())
    }
}
