use std::cmp::Ordering;

use crate::Cow;

impl PartialOrd for Cow<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
