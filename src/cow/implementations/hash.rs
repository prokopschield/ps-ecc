use std::hash::{Hash, Hasher};

use crate::Cow;

impl Hash for Cow<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (**self).hash(state);
    }
}

#[cfg(test)]
mod tests {
    use std::hash::{DefaultHasher, Hash, Hasher};

    use ps_buffer::Buffer;

    use crate::Cow;

    type TestError = Box<dyn std::error::Error>;

    fn hash(value: &Cow) -> u64 {
        let mut hasher = DefaultHasher::new();

        value.hash(&mut hasher);

        hasher.finish()
    }

    #[test]
    fn test_hash_agrees_across_variants() -> Result<(), TestError> {
        let borrowed = Cow::Borrowed(b"abc");
        let owned = Cow::Owned(Buffer::from_slice(b"abc")?.share());

        assert_eq!(hash(&borrowed), hash(&owned));

        Ok(())
    }
}
