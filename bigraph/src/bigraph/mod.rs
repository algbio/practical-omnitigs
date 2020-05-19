use std::ops::AddAssign;
use std::convert::{TryInto, TryFrom};
use std::fmt;

mod petgraph_impl;

pub trait StaticGraph<NodeData, EdgeData, IndexType> {
    type OutNeighbors: IntoIterator<Item = NodeIndex<IndexType>>;
    type InNeighbors: IntoIterator<Item = NodeIndex<IndexType>>;

    fn node_indices(&self) -> NodeIndices<IndexType>;

    fn edge_indices(&self) -> EdgeIndices<IndexType>;

    fn node_data(&self, node_id: NodeIndex<IndexType>) -> Option<&NodeData>;

    fn edge_data(&self, edge_id: EdgeIndex<IndexType>) -> Option<&EdgeData>;

    fn node_data_mut(&mut self, node_id: NodeIndex<IndexType>) -> Option<&mut NodeData>;

    fn edge_data_mut(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<&mut EdgeData>;

    // TODO This does not work as neighbors need a lifetime.
    fn out_neighbors(&self, node_id: NodeIndex<IndexType>) -> Option<Self::OutNeighbors>;

    fn in_neighbors(&self, node_id: NodeIndex<IndexType>) -> Option<Self::InNeighbors>;
}

pub trait DynamicGraph<NodeData, EdgeData, IndexType>: StaticGraph<NodeData, EdgeData, IndexType> {
    fn add_node(&mut self, node_data: NodeData) -> NodeIndex<IndexType>;

    fn add_edge(&mut self, edge_data: EdgeData) -> EdgeIndex<IndexType>;

    fn remove_node(&mut self, node_id: NodeIndex<IndexType>) -> Option<NodeData>;

    fn remove_edge(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<EdgeData>;
}

pub trait Bigraph<NodeData, EdgeData, IndexType>: StaticGraph<NodeData, EdgeData, IndexType> {
    fn reverse_complement_node(&self, node_id: NodeIndex<IndexType>) -> Option<NodeIndex<IndexType>>;

    fn reverse_complement_edge(&self, edge_id: EdgeIndex<IndexType>) -> Option<EdgeIndex<IndexType>>;
}

#[derive(Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy, OpaqueTypedef)]
#[opaque_typedef(derive(Display))]
pub struct NodeIndex<IndexType: Sized>(IndexType);
#[derive(Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy, OpaqueTypedef)]
#[opaque_typedef(derive(Display))]
pub struct EdgeIndex<IndexType: Sized>(IndexType);

pub struct NodeIndices<IndexType: Sized> {
    start: IndexType,
    end: IndexType,
}

pub struct EdgeIndices<IndexType: Sized> {
    start: IndexType,
    end: IndexType,
}

impl<RawType, IndexType: Sized + TryFrom<RawType> + fmt::Debug> From<(RawType, RawType)> for NodeIndices<IndexType>
    where <IndexType as std::convert::TryFrom<RawType>>::Error: std::fmt::Debug {
    fn from(raw: (RawType, RawType)) -> Self {
        Self {
            start: raw.0.try_into().unwrap(),
            end: raw.1.try_into().unwrap(),
        }
    }
}

impl<RawType, IndexType: Sized + TryFrom<RawType> + fmt::Debug> From<(RawType, RawType)> for EdgeIndices<IndexType>
    where <IndexType as std::convert::TryFrom<RawType>>::Error: std::fmt::Debug {
    fn from(raw: (RawType, RawType)) -> Self {
        Self {
            start: raw.0.try_into().unwrap(),
            end: raw.1.try_into().unwrap(),
        }
    }
}

impl<IndexType: Sized + Ord + Clone + AddAssign + From<u8>> Iterator for NodeIndices<IndexType> {
    type Item = NodeIndex<IndexType>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end {
            let result = Some(NodeIndex(self.start.clone()));
            self.start += 1.into();
            result
        } else {
            None
        }
    }
}

impl<IndexType: Sized + Ord + Clone + AddAssign + From<u8>> Iterator for EdgeIndices<IndexType> {
    type Item = EdgeIndex<IndexType>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end {
            let result = Some(EdgeIndex(self.start.clone()));
            self.start += 1.into();
            result
        } else {
            None
        }
    }
}