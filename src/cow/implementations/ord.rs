use std::cmp::Ordering;

use crate::Cow;

impl Ord for Cow<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        (**self).cmp(&**other)
    }
}

#[cfg(test)]
mod tests {
    use ps_buffer::Buffer;

    use crate::Cow;

    type TestError = Box<dyn std::error::Error>;

    #[test]
    fn test_ord_compares_content_across_variants() -> Result<(), TestError> {
        let borrowed = Cow::Borrowed(b"b");
        let owned = Cow::Owned(Buffer::from_slice(b"a")?.share());

        assert!(borrowed > owned);
        assert!(owned < borrowed);

        Ok(())
    }
}
