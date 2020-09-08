/// Traits concerning dynamic bigraphs, i.e. bigraphs that can be mutated.
pub mod dynamic_bigraph;
/// Traits concerning static bigraphs, i.e. bigraphs that are immutable.
pub mod static_bigraph;

/// A type of data that can be mirrored.
/// This is used to map nodes or edges to their mirrors.
pub trait BidirectedData {
    /// Compute the mirror of this data.
    fn mirror(&self) -> Self;
}

impl BidirectedData for () {
    fn mirror(&self) -> Self {}
}
