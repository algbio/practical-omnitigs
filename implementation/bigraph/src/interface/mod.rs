pub mod dynamic_bigraph;
pub mod static_bigraph;

pub trait BidirectedNodeData {
    fn reverse_complement(&self) -> Self;
}
