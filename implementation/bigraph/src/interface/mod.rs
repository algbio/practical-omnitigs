mod dynamic_bigraph;
mod static_bigraph;

pub use dynamic_bigraph::*;
pub use static_bigraph::*;

pub trait BidirectedNodeData {
    fn reverse_complement(&self) -> Self;
}
