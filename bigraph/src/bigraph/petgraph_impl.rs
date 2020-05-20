use petgraph::{Graph, Directed, Direction};
use crate::{ImmutableGraphContainer, EdgeIndex, NodeIndex, bigraph};
use crate::bigraph::{NodeIndices, EdgeIndices, NavigableGraph};
use petgraph::graph::Neighbors;
use std::iter::Map;

impl<NodeData, EdgeData> ImmutableGraphContainer<NodeData, EdgeData, usize> for Graph<NodeData, EdgeData, Directed, usize> {
    fn node_indices(&self) -> NodeIndices<usize> {
        NodeIndices::from((0, self.node_count()))
    }

    fn edge_indices(&self) -> EdgeIndices<usize> {
        EdgeIndices::from((0, self.edge_count()))
    }

    fn node_data(&self, node_id: NodeIndex<usize>) -> Option<&NodeData> {
        self.node_weight(node_id.0.into())
    }

    fn edge_data(&self, edge_id: EdgeIndex<usize>) -> Option<&EdgeData> {
        self.edge_weight(edge_id.0.into())
    }

    fn node_data_mut(&mut self, node_id: NodeIndex<usize>) -> Option<&mut NodeData> {
        self.node_weight_mut(node_id.0.into())
    }

    fn edge_data_mut(&mut self, edge_id: EdgeIndex<usize>) -> Option<&mut EdgeData> {
        self.edge_weight_mut(edge_id.0.into())
    }
}

impl<'a, NodeData, EdgeData> NavigableGraph<NodeData, EdgeData, usize> for &'a Graph<NodeData, EdgeData, Directed, usize> {
    type OutNeighbors = Map<Neighbors<'a, EdgeData, usize>, fn(petgraph::graph::NodeIndex<usize>) -> bigraph::NodeIndex<usize>>;
    //type OutNeighbors = Neighbors<'a, EdgeData, usize>;
    type InNeighbors = Map<Neighbors<'a, EdgeData, usize>, fn(petgraph::graph::NodeIndex<usize>) -> bigraph::NodeIndex<usize>>;

    fn out_neighbors(self, node_id: NodeIndex<usize>) -> Option<Self::OutNeighbors> {
        if node_id.0 < self.node_count() {
            Some(self.neighbors_directed(node_id.0.into(), Direction::Outgoing).map(|n| bigraph::NodeIndex(n.index())))
        } else {
            None
        }
    }

    fn in_neighbors(self, node_id: NodeIndex<usize>) -> Option<Self::InNeighbors> {
        if node_id.0 < self.node_count() {
            Some(self.neighbors_directed(node_id.0.into(), Direction::Incoming).map(|n| bigraph::NodeIndex(n.index())))
        } else {
            None
        }
    }
}

