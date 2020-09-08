use crate::index::GraphIndex;
use crate::interface::subgraph::Subgraph;
use crate::interface::ImmutableGraphContainer;
use bitvector::BitVector;

/// A subgraph that stores the presence or absence of a node or edge using bitvectors.
pub struct BitVectorSubgraph<'a, Graph> {
    original_graph: &'a Graph,
    present_nodes: BitVector,
    present_edges: BitVector,
}

impl<'a, Graph: ImmutableGraphContainer> Subgraph<'a, Graph> for BitVectorSubgraph<'a, Graph> {
    fn new_empty(graph: &'a Graph) -> Self {
        Self {
            original_graph: graph,
            present_nodes: BitVector::new(graph.node_count()),
            present_edges: BitVector::new(graph.edge_count()),
        }
    }

    fn new_full(graph: &'a Graph) -> Self {
        Self {
            original_graph: graph,
            present_nodes: BitVector::ones(graph.node_count()),
            present_edges: BitVector::ones(graph.edge_count()),
        }
    }

    fn original_graph(&self) -> &'a Graph {
        self.original_graph
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

    /// Returns the amount of nodes in the subgraph.
    fn node_count(&self) -> usize {
        self.present_nodes.len()
    }

    /// Returns the amount of edges in the subgraph.
    fn edge_count(&self) -> usize {
        self.present_edges.len()
    }
}