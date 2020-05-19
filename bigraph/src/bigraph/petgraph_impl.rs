use petgraph::{Graph, Directed};
use crate::{StaticGraph, EdgeIndex, NodeIndex};
use crate::bigraph::{NodeIndices, EdgeIndices};
use petgraph::graph::Neighbors;

impl<'a, NodeData, EdgeData> StaticGraph<NodeData, EdgeData, usize> for Graph<NodeData, EdgeData, Directed, usize> {
    type OutNeighbors = Neighbors<'a, EdgeData, usize>;
    type InNeighbors = ();

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

    fn out_neighbors(&self, node_id: NodeIndex<usize>) -> Option<Self::OutNeighbors> {
        unimplemented!()
    }

    fn in_neighbors(&self, node_id: NodeIndex<usize>) -> Option<Self::InNeighbors> {
        unimplemented!()
    }
}