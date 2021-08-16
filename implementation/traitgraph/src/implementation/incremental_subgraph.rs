use crate::index::GraphIndex;
use crate::interface::subgraph::DecoratingSubgraph;
use crate::interface::{GraphBase, ImmutableGraphContainer};

type IntegerType = usize;

/// A subgraph that stores the presence or absence of a node or edge using integers.
/// Additionally, this subgraph has a current step that can be altered.
/// Nodes and edges that are added are added with that step, and only nodes and edges with a step lower or equal to the current one are counted as present.
/// This allows to combined multiple subgraphs into one, if they are totally ordered by the subset relation.
pub struct IncrementalSubgraph<'a, Graph: GraphBase> {
    parent_graph: &'a Graph,
    present_nodes: Vec<IntegerType>,
    present_edges: Vec<IntegerType>,
    new_nodes: Vec<Vec<Graph::NodeIndex>>,
    new_edges: Vec<Vec<Graph::EdgeIndex>>,
    current_step: IntegerType,
}

impl<'a, Graph: ImmutableGraphContainer> IncrementalSubgraph<'a, Graph> {
    /// Create an incremental subgraph with the given amount of incremental steps.
    pub fn new_with_incremental_steps(
        graph: <Self as DecoratingSubgraph>::ParentGraphRef,
        incremental_steps: usize,
    ) -> Self {
        Self {
            parent_graph: graph,
            present_nodes: vec![IntegerType::max_value(); graph.node_count()],
            present_edges: vec![IntegerType::max_value(); graph.edge_count()],
            new_nodes: vec![Default::default(); incremental_steps],
            new_edges: vec![Default::default(); incremental_steps],
            current_step: 0,
        }
    }

    /// Set the current incremental step of the graph.
    pub fn set_current_step(&mut self, current_step: IntegerType) {
        debug_assert!(current_step < self.new_nodes.len() && current_step < self.new_edges.len());
        self.current_step = current_step;
    }

    /// Return the nodes that are added in the current incremental step.
    pub fn new_nodes(&self) -> &Vec<Graph::NodeIndex> {
        debug_assert!(self.current_step < self.new_nodes.len());
        &self.new_nodes[self.current_step]
    }

    /// Return the edges that are added in the current incremental step.
    pub fn new_edges(&self) -> &Vec<Graph::EdgeIndex> {
        debug_assert!(self.current_step < self.new_edges.len());
        &self.new_edges[self.current_step]
    }
}

impl<'a, Graph: ImmutableGraphContainer> DecoratingSubgraph for IncrementalSubgraph<'a, Graph> {
    type ParentGraph = Graph;
    type ParentGraphRef = &'a Graph;

    /// Not implemented for this type.
    fn new_empty(_graph: Self::ParentGraphRef) -> Self {
        /*Self {
            parent_graph: graph,
            present_nodes: vec![0; graph.node_count()],
            present_edges: vec![0; graph.edge_count()],
            new_nodes: Default::default(),
            new_edges: Default::default(),
            current_index: 0,
        }*/
        unimplemented!()
    }

    /// Not implemented for this type.
    fn new_full(_graph: Self::ParentGraphRef) -> Self {
        unimplemented!()
    }

    fn clear(&mut self) {
        for node in &mut self.present_nodes {
            *node = IntegerType::max_value();
        }
        for edge in &mut self.present_edges {
            *edge = IntegerType::max_value();
        }
        for nodes in &mut self.new_nodes {
            nodes.clear();
        }
        for edges in &mut self.new_edges {
            edges.clear();
        }
    }

    fn parent_graph(&self) -> &Self::ParentGraph {
        self.parent_graph
    }

    fn contains_node(&self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex) -> bool {
        debug_assert!(node_index.as_usize() < self.present_nodes.capacity());
        self.present_nodes[node_index.as_usize()] <= self.current_step
    }

    fn contains_edge(&self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex) -> bool {
        debug_assert!(edge_index.as_usize() < self.present_edges.capacity());
        self.present_edges[edge_index.as_usize()] <= self.current_step
    }

    /// Panics if the node_index is not valid for the graph passed in the constructor.
    /// Panics also if the node was added already.
    fn add_node(&mut self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex) {
        debug_assert!(node_index.as_usize() < self.present_nodes.capacity());
        if self.present_nodes[node_index.as_usize()] > self.current_step {
            debug_assert_eq!(
                self.present_nodes[node_index.as_usize()],
                IntegerType::max_value()
            );
            self.present_nodes[node_index.as_usize()] = self.current_step;
            self.new_nodes[self.current_step].push(node_index);
        }
    }

    /// Panics if the edge_index is not valid for the graph passed in the constructor.
    /// Panics also if the edge was added already.
    fn add_edge(&mut self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex) {
        debug_assert!(edge_index.as_usize() < self.present_edges.capacity());
        if self.present_edges[edge_index.as_usize()] > self.current_step {
            debug_assert_eq!(
                self.present_edges[edge_index.as_usize()],
                IntegerType::max_value()
            );
            self.present_edges[edge_index.as_usize()] = self.current_step;
            self.new_edges[self.current_step].push(edge_index);
        }
    }

    /// Panics if the node_index is not valid for the graph passed in the constructor.
    fn remove_node(&mut self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex) {
        debug_assert!(node_index.as_usize() < self.present_nodes.capacity());
        self.present_nodes[node_index.as_usize()] =
            self.present_nodes[node_index.as_usize()].max(self.current_step + 1);
        unimplemented!("This method is not intended to be called on an incremental subgraph.");
    }

    /// Panics if the edge_index is not valid for the graph passed in the constructor.
    fn remove_edge(&mut self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex) {
        debug_assert!(edge_index.as_usize() < self.present_edges.capacity());
        self.present_edges[edge_index.as_usize()] =
            self.present_edges[edge_index.as_usize()].max(self.current_step + 1);
        unimplemented!("This method is not intended to be called on an incremental subgraph.");
    }

    /// Returns the amount of nodes in the subgraph.
    fn node_count(&self) -> usize {
        self.present_nodes
            .iter()
            .filter(|&&n| n <= self.current_step)
            .count()
    }

    /// Returns the amount of edges in the subgraph.
    fn edge_count(&self) -> usize {
        self.present_edges
            .iter()
            .filter(|&&n| n <= self.current_step)
            .count()
    }
}
