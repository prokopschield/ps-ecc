mod implementations;
mod methods;

/// number of distinct codewords covering each byte
#[derive(Clone, Copy, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum OverlapFactor {
    /// each byte is covered once (plain Reed-Solomon, no overlap)
    #[default]
    Simple = 0,

    /// each byte is covered twice
    Double = 64,

    /// each byte is covered three times
    Triple = 128,

    /// each byte is covered four times
    Quadruple = 192,
}
