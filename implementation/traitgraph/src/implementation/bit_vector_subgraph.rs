use crate::index::GraphIndex;
use crate::interface::subgraph::Subgraph;
use crate::interface::ImmutableGraphContainer;
use bitvector::BitVector;

/// A subgraph that stores the presence or absence of a node or edge using bitvectors.
pub struct BitVectorSubgraph {
    present_nodes: BitVector,
    present_edges: BitVector,
}

impl<Graph: ImmutableGraphContainer> Subgraph<Graph> for BitVectorSubgraph {
    fn new_empty(graph: &Graph) -> Self {
        Self {
            present_nodes: BitVector::new(graph.node_count()),
            present_edges: BitVector::new(graph.edge_count()),
        }
    }

    fn new_full(graph: &Graph) -> Self {
        Self {
            present_nodes: BitVector::ones(graph.node_count()),
            present_edges: BitVector::ones(graph.edge_count()),
        }
    }

    fn contains_node(&self, node_index: Graph::NodeIndex) -> bool {
        assert!(node_index.as_usize() < self.present_nodes.capacity());
        self.present_nodes.contains(node_index.as_usize())
    }

    fn contains_edge(&self, edge_index: Graph::EdgeIndex) -> bool {
        assert!(edge_index.as_usize() < self.present_edges.capacity());
        self.present_edges.contains(edge_index.as_usize())
    }

    /// Panics if the node_index is not valid for the graph passed in the constructor.
    fn add_node(&mut self, node_index: Graph::NodeIndex) {
        assert!(node_index.as_usize() < self.present_nodes.capacity());
        self.present_nodes.insert(node_index.as_usize());
    }

    /// Panics if the edge_index is not valid for the graph passed in the constructor.
    fn add_edge(&mut self, edge_index: Graph::EdgeIndex) {
        assert!(edge_index.as_usize() < self.present_edges.capacity());
        self.present_edges.insert(edge_index.as_usize());
    }

    /// Panics if the node_index is not valid for the graph passed in the constructor.
    fn remove_node(&mut self, node_index: Graph::NodeIndex) {
        assert!(node_index.as_usize() < self.present_nodes.capacity());
        self.present_nodes.remove(node_index.as_usize());
    }

    /// Panics if the edge_index is not valid for the graph passed in the constructor.
    fn remove_edge(&mut self, edge_index: Graph::EdgeIndex) {
        assert!(edge_index.as_usize() < self.present_edges.capacity());
        self.present_edges.remove(edge_index.as_usize());
    }
}
