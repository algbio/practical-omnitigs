use crate::index::GraphIndex;
use crate::interface::subgraph::DecoratingSubgraph;
use crate::interface::{GraphBase, ImmutableGraphContainer};
use bitvector::BitVector;

/// A subgraph that stores the presence or absence of a node or edge using bitvectors.
pub struct BitVectorSubgraph<'a, Graph> {
    parent_graph: &'a Graph,
    present_nodes: BitVector,
    present_edges: BitVector,
}

impl<'a, Graph: ImmutableGraphContainer> DecoratingSubgraph for BitVectorSubgraph<'a, Graph> {
    type ParentGraph = Graph;
    type ParentGraphRef = &'a Graph;

    fn new_empty(graph: Self::ParentGraphRef) -> Self {
        Self {
            parent_graph: graph,
            present_nodes: BitVector::new(graph.node_count()),
            present_edges: BitVector::new(graph.edge_count()),
        }
    }

    fn new_full(graph: Self::ParentGraphRef) -> Self {
        Self {
            parent_graph: graph,
            present_nodes: BitVector::ones(graph.node_count()),
            present_edges: BitVector::ones(graph.edge_count()),
        }
    }

    fn clear(&mut self) {
        self.present_nodes.clear();
        self.present_edges.clear();
    }

    fn parent_graph(&self) -> &Self::ParentGraph {
        self.parent_graph
    }

    fn contains_node(&self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex) -> bool {
        debug_assert!(node_index.as_usize() < self.present_nodes.capacity());
        self.present_nodes.contains(node_index.as_usize())
    }

    fn contains_edge(&self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex) -> bool {
        debug_assert!(edge_index.as_usize() < self.present_edges.capacity());
        self.present_edges.contains(edge_index.as_usize())
    }

    /// Panics if the node_index is not valid for the graph passed in the constructor.
    fn add_node(&mut self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex) {
        debug_assert!(node_index.as_usize() < self.present_nodes.capacity());
        self.present_nodes.insert(node_index.as_usize());
    }

    /// Panics if the edge_index is not valid for the graph passed in the constructor.
    fn add_edge(&mut self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex) {
        debug_assert!(edge_index.as_usize() < self.present_edges.capacity());
        self.present_edges.insert(edge_index.as_usize());
    }

    /// Panics if the node_index is not valid for the graph passed in the constructor.
    fn remove_node(&mut self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex) {
        debug_assert!(node_index.as_usize() < self.present_nodes.capacity());
        self.present_nodes.remove(node_index.as_usize());
    }

    /// Panics if the edge_index is not valid for the graph passed in the constructor.
    fn remove_edge(&mut self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex) {
        debug_assert!(edge_index.as_usize() < self.present_edges.capacity());
        self.present_edges.remove(edge_index.as_usize());
    }

    /// Returns the amount of nodes in the subgraph.
    fn node_count(&self) -> usize {
        self.present_nodes.len()
    }

    /// Returns the amount of edges in the subgraph.
    fn edge_count(&self) -> usize {
        self.present_edges.len()
    }
}
