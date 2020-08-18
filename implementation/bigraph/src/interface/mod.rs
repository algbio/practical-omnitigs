pub mod dynamic_bigraph;
pub mod static_bigraph;

pub trait BidirectedData {
    fn reverse_complement(&self) -> Self;
}

impl BidirectedData for () {
    fn reverse_complement(&self) -> Self {}
}
