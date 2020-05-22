use petgraph::{Graph, Directed, Direction};
use crate::{ImmutableGraphContainer, EdgeIndex, NodeIndex};
use crate::{NodeIndices, EdgeIndices, NavigableGraph};
use petgraph::graph::Neighbors;
use std::iter::Map;
use num_traits::{PrimInt, ToPrimitive};

impl<NodeData, EdgeData> ImmutableGraphContainer<NodeData, EdgeData, usize> for Graph<NodeData, EdgeData, Directed, usize> {
    fn node_indices(&self) -> NodeIndices<usize> {
        NodeIndices::from((0, self.node_count()))
    }

    fn edge_indices(&self) -> EdgeIndices<usize> {
        EdgeIndices::from((0, self.edge_count()))
    }

    fn node_count(&self) -> usize {
        self.node_count()
    }

    fn edge_count(&self) -> usize {
        self.edge_count()
    }

    fn node_data(&self, node_id: NodeIndex<usize>) -> Option<&NodeData> {
        self.node_weight(node_id.into())
    }

    fn edge_data(&self, edge_id: EdgeIndex<usize>) -> Option<&EdgeData> {
        self.edge_weight(edge_id.into())
    }

    fn node_data_mut(&mut self, node_id: NodeIndex<usize>) -> Option<&mut NodeData> {
        self.node_weight_mut(node_id.into())
    }

    fn edge_data_mut(&mut self, edge_id: EdgeIndex<usize>) -> Option<&mut EdgeData> {
        self.edge_weight_mut(edge_id.into())
    }
}

impl<'a, NodeData, EdgeData: 'a> NavigableGraph<'a, NodeData, EdgeData, usize> for Graph<NodeData, EdgeData, Directed, usize> {
    type OutNeighbors = Map<Neighbors<'a, EdgeData, usize>, fn(petgraph::graph::NodeIndex<usize>) -> NodeIndex<usize>>;
    //type OutNeighbors = Neighbors<'a, EdgeData, usize>;
    type InNeighbors = Map<Neighbors<'a, EdgeData, usize>, fn(petgraph::graph::NodeIndex<usize>) -> NodeIndex<usize>>;

    fn out_neighbors(&'a self, node_id: NodeIndex<usize>) -> Option<Self::OutNeighbors> {
        if node_id < self.node_count().into() {
            Some(self.neighbors_directed(node_id.into(), Direction::Outgoing).map(|n| NodeIndex::from(n.index())))
        } else {
            None
        }
    }

    fn in_neighbors(&'a self, node_id: NodeIndex<usize>) -> Option<Self::InNeighbors> {
        if node_id < self.node_count().into() {
            Some(self.neighbors_directed(node_id.into(), Direction::Incoming).map(|n| NodeIndex::from(n.index())))
        } else {
            None
        }
    }
}

impl<IndexType: PrimInt + ToPrimitive + petgraph::graph::IndexType> From<NodeIndex<IndexType>> for petgraph::graph::NodeIndex<IndexType> {
    fn from(index: NodeIndex<IndexType>) -> Self {
        IndexType::from(index).unwrap().into()
    }
}

impl<IndexType: PrimInt + ToPrimitive + petgraph::graph::IndexType> From<EdgeIndex<IndexType>> for petgraph::graph::EdgeIndex<IndexType> {
    fn from(index: EdgeIndex<IndexType>) -> Self {
        IndexType::from(index).unwrap().into()
    }
}