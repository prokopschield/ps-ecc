mod implementations;
mod methods;

/// Approximate number of distinct codewords covering each byte.
///
/// Integer division of the segment stride makes coverage uneven for the
/// overlapping factors: bytes near segment boundaries may be covered one
/// time more, and leading bytes fewer times, than the nominal count.
#[derive(Clone, Copy, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum OverlapFactor {
    /// Each byte is covered once (plain Reed-Solomon, no overlap).
    #[default]
    Simple = 0,

    /// Each byte is covered approximately twice.
    Double = 64,

    /// Each byte is covered approximately three times.
    Triple = 128,

    /// Each byte is covered approximately four times.
    Quadruple = 192,
}
