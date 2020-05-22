use std::ops::AddAssign;
use std::convert::{TryInto, TryFrom};
use std::fmt;
use num_traits::{PrimInt, ToPrimitive, NumCast};
use crate::{NodeIndices, EdgeIndices, NodeIndex, EdgeIndex};

pub trait ImmutableGraphContainer<NodeData, EdgeData, IndexType: PrimInt> {
    fn node_indices(&self) -> NodeIndices<IndexType>;

    fn edge_indices(&self) -> EdgeIndices<IndexType>;

    fn node_count(&self) -> IndexType;

    fn edge_count(&self) -> IndexType;

    fn node_data(&self, node_id: NodeIndex<IndexType>) -> Option<&NodeData>;

    fn edge_data(&self, edge_id: EdgeIndex<IndexType>) -> Option<&EdgeData>;

    fn node_data_mut(&mut self, node_id: NodeIndex<IndexType>) -> Option<&mut NodeData>;

    fn edge_data_mut(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<&mut EdgeData>;
}

pub trait MutableGraphContainer<NodeData, EdgeData, IndexType> {
    fn add_node(&mut self, node_data: NodeData) -> NodeIndex<IndexType>;

    fn add_edge(&mut self, edge_data: EdgeData) -> EdgeIndex<IndexType>;

    fn remove_node(&mut self, node_id: NodeIndex<IndexType>) -> Option<NodeData>;

    fn remove_edge(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<EdgeData>;
}

pub trait NavigableGraph<'a, NodeData, EdgeData, IndexType> {
    type OutNeighbors: IntoIterator<Item = NodeIndex<IndexType>>;
    type InNeighbors: IntoIterator<Item = NodeIndex<IndexType>>;

    fn out_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::OutNeighbors>;

    fn in_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::InNeighbors>;
}

pub trait Bigraph<NodeData, EdgeData, IndexType> {
    fn reverse_complement_node(&self, node_id: NodeIndex<IndexType>) -> Option<NodeIndex<IndexType>>;

    fn reverse_complement_edge(&self, edge_id: EdgeIndex<IndexType>) -> Option<EdgeIndex<IndexType>>;
}

pub trait StaticGraph<NodeData, EdgeData, IndexType: PrimInt>: ImmutableGraphContainer<NodeData, EdgeData, IndexType> + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType> {}
impl<NodeData, EdgeData, IndexType: PrimInt, T: ImmutableGraphContainer<NodeData, EdgeData, IndexType> + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType>> StaticGraph<NodeData, EdgeData, IndexType> for T {}

pub trait DynamicGraph<NodeData, EdgeData, IndexType: PrimInt>: ImmutableGraphContainer<NodeData, EdgeData, IndexType> + MutableGraphContainer<NodeData, EdgeData, IndexType> + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType> {}
impl<NodeData, EdgeData, IndexType: PrimInt, T: ImmutableGraphContainer<NodeData, EdgeData, IndexType> + MutableGraphContainer<NodeData, EdgeData, IndexType> + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType>> DynamicGraph<NodeData, EdgeData, IndexType> for T {}

pub trait StaticBigraph<NodeData, EdgeData, IndexType: PrimInt>: ImmutableGraphContainer<NodeData, EdgeData, IndexType> + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType> + Bigraph<NodeData, EdgeData, IndexType> {}
impl<NodeData, EdgeData, IndexType: PrimInt, T: ImmutableGraphContainer<NodeData, EdgeData, IndexType> + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType> + Bigraph<NodeData, EdgeData, IndexType>> StaticBigraph<NodeData, EdgeData, IndexType> for T {}
