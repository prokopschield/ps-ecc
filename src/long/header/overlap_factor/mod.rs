mod implementations;
mod methods;

/// approximate number of distinct codewords covering each byte
///
/// Integer division of the segment stride makes coverage uneven for the
/// overlapping factors: bytes near segment boundaries may be covered one
/// time more, and leading bytes fewer times, than the nominal count.
#[derive(Clone, Copy, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum OverlapFactor {
    /// each byte is covered once (plain Reed-Solomon, no overlap)
    #[default]
    Simple = 0,

    /// each byte is covered approximately twice
    Double = 64,

    /// each byte is covered approximately three times
    Triple = 128,

    /// each byte is covered approximately four times
    Quadruple = 192,
}
