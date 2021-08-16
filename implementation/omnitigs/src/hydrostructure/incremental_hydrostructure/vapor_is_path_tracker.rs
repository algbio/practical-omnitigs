use crate::hydrostructure::incremental_hydrostructure::IncrementalSafetyTracker;
use traitgraph::implementation::bit_vector_subgraph::BitVectorSubgraph;
use traitgraph::implementation::incremental_subgraph::IncrementalSubgraph;
use traitgraph::index::GraphIndex;
use traitgraph::interface::subgraph::DecoratingSubgraph;
use traitgraph::interface::ImmutableGraphContainer;
use traitgraph::interface::{GraphBase, NavigableGraph, StaticGraph};

/// A structure to dynamically track if the vapor is a path.
/// This structure takes O(|Edges|) time to construct and offers an O(1) query if the subgraph is a path.
/// Inserting and removing each node and edge of the graph once in any order takes O(|Edges|) time.
///
/// All runtimes assume that the endpoints of an edge can be retrieved in O(1).
pub struct VaporIsPathTracker<'a, Graph> {
    subgraph: BitVectorSubgraph<'a, Graph>,
    node_in_edges: Vec<usize>,
    node_out_edges: Vec<usize>,
    source_count: usize,
    sink_count: usize,
    split_node_count: usize,
    join_node_count: usize,
}

impl<'a, Graph: StaticGraph> IncrementalSafetyTracker<'a, Graph> for VaporIsPathTracker<'a, Graph> {
    fn new_with_empty_subgraph(graph: &'a Graph) -> Self {
        Self {
            node_in_edges: vec![0; graph.node_count()],
            node_out_edges: vec![0; graph.node_count()],
            subgraph: BitVectorSubgraph::new_empty(graph),
            source_count: 0,
            sink_count: 0,
            split_node_count: 0,
            join_node_count: 0,
        }
    }

    fn clear(&mut self) {
        if DecoratingSubgraph::node_count(&self.subgraph) == 0
            && DecoratingSubgraph::edge_count(&self.subgraph) == 0
        {
            return;
        }

        self.subgraph.clear();
        for in_edge in &mut self.node_in_edges {
            *in_edge = 0;
        }
        for out_edge in &mut self.node_out_edges {
            *out_edge = 0;
        }
        self.source_count = 0;
        self.sink_count = 0;
        self.split_node_count = 0;
        self.join_node_count = 0;
    }

    fn reset(&mut self, r_plus: &IncrementalSubgraph<Graph>, r_minus: &IncrementalSubgraph<Graph>) {
        self.clear();
        for node in r_plus.parent_graph().node_indices() {
            if r_plus.contains_node(node) && r_minus.contains_node(node) {
                self.add_node(node);
            }
        }
        for edge in r_plus.parent_graph().edge_indices() {
            if r_plus.contains_edge(edge) && r_minus.contains_edge(edge) {
                self.add_edge(edge);
            }
        }
    }

    fn add_incremental_subgraph_step(
        &mut self,
        r_plus: &IncrementalSubgraph<Graph>,
        r_minus: &IncrementalSubgraph<Graph>,
    ) {
        for node in r_plus.new_nodes() {
            if r_minus.contains_node(*node) {
                self.add_node(*node);
            }
        }
        for edge in r_plus.new_edges() {
            if r_minus.contains_edge(*edge) {
                self.add_edge(*edge);
            }
        }
    }

    fn remove_incremental_subgraph_step(
        &mut self,
        _r_plus: &IncrementalSubgraph<Graph>,
        r_minus: &IncrementalSubgraph<Graph>,
    ) {
        for node in r_minus.new_nodes() {
            if self.contains_node(*node) {
                self.remove_node(*node);
            }
        }
        for edge in r_minus.new_edges() {
            if self.contains_edge(*edge) {
                self.remove_edge(*edge);
            }
        }
    }

    /// Returns true if the subgraph is a path.
    /// Empty graphs are not counted as paths.
    fn is_safe(&self, is_forward_univocal: bool, is_backward_univocal: bool) -> bool {
        is_forward_univocal
            || is_backward_univocal
            || (self.source_count == 1
                && self.sink_count == 1
                && self.split_node_count == 0
                && self.join_node_count == 0)
    }

    fn does_safety_equal_bridge_like() -> bool {
        true
    }
}

impl<'a, Graph: StaticGraph> VaporIsPathTracker<'a, Graph> {
    /// Returns true if the given node is in the subgraph.
    pub fn contains_node(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        self.subgraph.contains_node(node)
    }

    /// Returns true if the given edge is in the subgraph.
    pub fn contains_edge(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        DecoratingSubgraph::contains_edge(&self.subgraph, edge)
    }

    /// Add a node to the subgraph.
    pub fn add_node(&mut self, node: <Graph as GraphBase>::NodeIndex) {
        debug_assert!(!self.subgraph.contains_node(node));
        self.subgraph.add_node(node);

        let out_degree = self.subgraph.out_degree(node);
        let in_degree = self.subgraph.in_degree(node);

        self.node_in_edges[node.as_usize()] = in_degree;
        self.node_out_edges[node.as_usize()] = out_degree;

        // Remove existing out (in) edges from sources (sinks).
        self.source_count -= out_degree;
        self.sink_count -= in_degree;

        // If the new node is a source or sink, update accordingly.
        if in_degree == 0 {
            self.source_count += 1;
        }
        if out_degree == 0 {
            self.sink_count += 1;
        }

        // If the new node is a split or join, update accordingly.
        if in_degree > 1 {
            self.join_node_count += 1;
        }
        if out_degree > 1 {
            self.split_node_count += 1;
        }
    }

    /// Add an edge to the subgraph.
    pub fn add_edge(&mut self, edge: <Graph as GraphBase>::EdgeIndex) {
        debug_assert!(
            !DecoratingSubgraph::contains_edge(&self.subgraph, edge),
            "Subgraph already contains edge {:?}",
            edge
        );
        self.subgraph.add_edge(edge);

        let endpoints = self.subgraph.edge_endpoints(edge);
        let from_node = endpoints.from_node;
        let to_node = endpoints.to_node;

        if self.subgraph.contains_node(from_node) {
            self.node_out_edges[from_node.as_usize()] += 1;

            match self.node_out_edges[from_node.as_usize()] {
                1 => self.sink_count -= 1,
                2 => self.split_node_count += 1,
                _ => {}
            }
        } else {
            self.source_count += 1;
        }

        if self.subgraph.contains_node(to_node) {
            self.node_in_edges[to_node.as_usize()] += 1;

            match self.node_in_edges[to_node.as_usize()] {
                1 => self.source_count -= 1,
                2 => self.join_node_count += 1,
                _ => {}
            }
        } else {
            self.sink_count += 1;
        }
    }

    /// Remove a node from the subgraph.
    pub fn remove_node(&mut self, node: <Graph as GraphBase>::NodeIndex) {
        debug_assert!(self.subgraph.contains_node(node));

        let out_degree = self.subgraph.out_degree(node);
        let in_degree = self.subgraph.in_degree(node);
        self.subgraph.remove_node(node);

        self.node_in_edges[node.as_usize()] = 0;
        self.node_out_edges[node.as_usize()] = 0;

        // Add existing out (in) edges to sources (sinks).
        self.source_count += out_degree;
        self.sink_count += in_degree;

        // If the deleted node is a source or sink, update accordingly.
        if in_degree == 0 {
            self.source_count -= 1;
        }
        if out_degree == 0 {
            self.sink_count -= 1;
        }

        // If the deleted node is a split or join, update accordingly.
        if in_degree > 1 {
            self.join_node_count -= 1;
        }
        if out_degree > 1 {
            self.split_node_count -= 1;
        }
    }

    /// Remove an edge from the subgraph.
    pub fn remove_edge(&mut self, edge: <Graph as GraphBase>::EdgeIndex) {
        debug_assert!(DecoratingSubgraph::contains_edge(&self.subgraph, edge));

        let endpoints = self.subgraph.edge_endpoints(edge);
        self.subgraph.remove_edge(edge);
        let from_node = endpoints.from_node;
        let to_node = endpoints.to_node;

        if self.subgraph.contains_node(from_node) {
            self.node_out_edges[from_node.as_usize()] -= 1;

            match self.node_out_edges[from_node.as_usize()] {
                0 => self.sink_count += 1,
                1 => self.split_node_count -= 1,
                _ => {}
            }
        } else {
            self.source_count -= 1;
        }

        if self.subgraph.contains_node(to_node) {
            self.node_in_edges[to_node.as_usize()] -= 1;

            match self.node_in_edges[to_node.as_usize()] {
                0 => self.source_count += 1,
                1 => self.join_node_count -= 1,
                _ => {}
            }
        } else {
            self.sink_count -= 1;
        }
    }
}

impl<'a, Graph: ImmutableGraphContainer> std::fmt::Debug for VaporIsPathTracker<'a, Graph>
where
    Graph::NodeIndex: std::fmt::Debug,
    Graph::EdgeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "VaporIsPathTracker[nodes: [")?;

        let mut once = true;
        for node in self.subgraph.parent_graph().node_indices() {
            if self.subgraph.contains_node(node) {
                if once {
                    once = false;
                } else {
                    write!(f, ", ")?;
                }
                write!(f, "{:?}", node)?;
            }
        }

        write!(f, "], edges: [")?;

        let mut once = true;
        for edge in self.subgraph.parent_graph().edge_indices() {
            if self.subgraph.contains_edge(edge) {
                if once {
                    once = false;
                } else {
                    write!(f, ", ")?;
                }
                write!(f, "{:?}", edge)?;
            }
        }

        write!(f, "]]")
    }
}
