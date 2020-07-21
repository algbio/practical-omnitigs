use crate::{EdgeIndex, EdgeIndices, NodeIndex, NodeIndices};
use num_traits::PrimInt;

mod dynamic_bigraph;
mod static_bigraph;

pub use dynamic_bigraph::*;
pub use static_bigraph::*;

pub struct Edge<'a, EdgeData, IndexType: PrimInt> {
    pub from_node: NodeIndex<IndexType>,
    pub to_node: NodeIndex<IndexType>,
    pub edge_data: &'a EdgeData,
}

pub trait ImmutableGraphContainer<'a, NodeData, EdgeData: 'a, IndexType: PrimInt> {
    type EdgeIterator: IntoIterator<Item = Edge<'a, EdgeData, IndexType>>;

    fn node_indices(&self) -> NodeIndices<IndexType>;

    fn edge_indices(&self) -> EdgeIndices<IndexType>;

    fn edges(&'a self) -> Self::EdgeIterator;

    fn node_count(&self) -> usize;

    fn edge_count(&self) -> usize;

    fn node_data(&self, node_id: NodeIndex<IndexType>) -> Option<&NodeData>;

    fn edge_data(&self, edge_id: EdgeIndex<IndexType>) -> Option<&EdgeData>;

    fn node_data_mut(&mut self, node_id: NodeIndex<IndexType>) -> Option<&mut NodeData>;

    fn edge_data_mut(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<&mut EdgeData>;

    fn contains_edge(&self, from: NodeIndex<IndexType>, to: NodeIndex<IndexType>) -> bool;
}

pub trait MutableGraphContainer<NodeData, EdgeData, IndexType> {
    fn add_node(&mut self, node_data: NodeData) -> NodeIndex<IndexType>;

    fn add_edge(
        &mut self,
        from: NodeIndex<IndexType>,
        to: NodeIndex<IndexType>,
        edge_data: EdgeData,
    ) -> EdgeIndex<IndexType>;

    fn remove_node(&mut self, node_id: NodeIndex<IndexType>) -> Option<NodeData>;

    fn remove_edge(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<EdgeData>;
}

pub trait NavigableGraph<'a, NodeData, EdgeData, IndexType> {
    type OutNeighbors: IntoIterator<Item = Neighbor<IndexType>>;
    type InNeighbors: IntoIterator<Item = Neighbor<IndexType>>;

    fn out_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::OutNeighbors>;

    fn in_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::InNeighbors>;
}

pub trait BidirectedNodeData {
    fn reverse_complement(&self) -> Self;
}

pub trait StaticGraph<'a, NodeData, EdgeData: 'a, IndexType: PrimInt>:
    ImmutableGraphContainer<'a, NodeData, EdgeData, IndexType>
    + for<'b> NavigableGraph<'b, NodeData, EdgeData, IndexType>
{
}
impl<
        'a,
        NodeData,
        EdgeData: 'a,
        IndexType: PrimInt,
        T: ImmutableGraphContainer<'a, NodeData, EdgeData, IndexType>
            + for<'b> NavigableGraph<'b, NodeData, EdgeData, IndexType>,
    > StaticGraph<'a, NodeData, EdgeData, IndexType> for T
{
}

pub trait DynamicGraph<'a, NodeData, EdgeData: 'a, IndexType: PrimInt>:
    StaticGraph<'a, NodeData, EdgeData, IndexType>
    + MutableGraphContainer<NodeData, EdgeData, IndexType>
{
}
impl<
        'a,
        NodeData,
        EdgeData: 'a,
        IndexType: PrimInt,
        T: StaticGraph<'a, NodeData, EdgeData, IndexType>
            + MutableGraphContainer<NodeData, EdgeData, IndexType>,
    > DynamicGraph<'a, NodeData, EdgeData, IndexType> for T
{
}

/*pub struct Edge<IndexType> {
    pub from: NodeIndex<IndexType>,
    pub to: NodeIndex<IndexType>,
}*/

pub struct Neighbor<IndexType> {
    pub edge_id: EdgeIndex<IndexType>,
    pub node_id: NodeIndex<IndexType>,
}
