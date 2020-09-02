use crate::interface::ImmutableGraphContainer;

/// A subset of nodes and edges of a graph.
/// This trait operates only on indices without being coupled to a graph instance.
/// There is not restriction on which nodes or edges must be contained in combination, so e.g. an edge from node 1 to node 4 can be contained, even if node 1 and node 4 are both not contained.
pub trait Subgraph<Graph: ImmutableGraphContainer> {
    /// Constructs a subgraph from the given graph without any nodes or edges.
    /// If not defined otherwise in the implementation, all node and edge ids of the given graph are valid arguments for the methods of this trait on the returned object.
    fn new_empty(graph: &Graph) -> Self;

    /// Constructs a subgraph from the given graph with all nodes and edges.
    /// If not defined otherwise in the implementation, all node and edge ids of the given graph are valid arguments for the methods of this trait on the returned object.
    fn new_full(graph: &Graph) -> Self;

    /// Returns true if the given node id is part of the subgraph.
    fn contains_node(&self, node_index: Graph::NodeIndex) -> bool;

    /// Returns true if the given edge id is part of the subgraph.
    fn contains_edge(&self, edge_index: Graph::EdgeIndex) -> bool;

    /// Adds the given node id to the subgraph.
    /// Note that some implementations may require an initial set of nodes to be known when a subgraph is created, and inserting nodes outside of this initial set might panic.
    fn add_node(&mut self, node_index: Graph::NodeIndex);

    /// Adds the given edge id to the subgraph.
    /// Note that some implementations may require an initial set of nodes to be known when a subgraph is created, and inserting nodes outside of this initial set might panic.
    fn add_edge(&mut self, edge_index: Graph::EdgeIndex);

    /// Removes the given node id from the graph.
    fn remove_node(&mut self, node_index: Graph::NodeIndex);

    /// Removes the given edge id from the graph.
    fn remove_edge(&mut self, edge_index: Graph::EdgeIndex);
}
