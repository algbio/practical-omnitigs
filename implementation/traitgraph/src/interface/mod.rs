//! The graph traits.
//!
//! The traits are roughly split up by different access types:
//!  - immutable reference (`ImmutableGraphContainer`)
//!  - mutable reference (`MutableGraphContainer`)
//!  - immutable reference that must outlive the return value (`TraversableGraph`)
//!
//! As it happens, the access types match well to common graph use cases, i.e. queries for nodes and edges, adding and removing nodes and edges as well as iterating over the neighbors of a node.

use crate::index::{GraphIndex, GraphIndices, OptionalGraphIndex};
use crate::walks::{EdgeWalk, NodeWalk};
use std::iter::FromIterator;

/// A set of traits for subgraphs.
/// A subgraph is a graph that is backed by an actual graph implementation, but that filters out some nodes or edges.
pub mod subgraph;

/// Contains the associated types of a graph.
pub trait GraphBase {
    /// The data type associated with each node.
    type NodeData;
    /// The data type associated with each edge.
    type EdgeData;
    /// The optional index type used for nodes.
    type OptionalNodeIndex: OptionalGraphIndex<Self::NodeIndex>;
    /// The optional index type used for edges.
    type OptionalEdgeIndex: OptionalGraphIndex<Self::EdgeIndex>;
    /// The index type used for nodes.
    type NodeIndex: GraphIndex<Self::OptionalNodeIndex>;
    /// The index type used for edges.
    type EdgeIndex: GraphIndex<Self::OptionalEdgeIndex>;

    /// Returns the none value of the optional node index type used by the trait.
    fn new_none_optional_node_index(&self) -> Self::OptionalNodeIndex {
        Self::OptionalNodeIndex::new_none()
    }

    /// Returns the none value of the optional edge index type used by the trait.
    fn new_none_optional_edge_index(&self) -> Self::OptionalEdgeIndex {
        Self::OptionalEdgeIndex::new_none()
    }
}

/// A container that contains a set of nodes and edges.
///
/// Graphs that implement this trait must have their nodes and edges indexed consecutively.
pub trait ImmutableGraphContainer: GraphBase {
    /// Returns an iterator over the node indices in this graph.
    fn node_indices(&self) -> GraphIndices<Self::NodeIndex, Self::OptionalNodeIndex>;

    /// Returns an iterator over the edge indices in this graph.
    fn edge_indices(&self) -> GraphIndices<Self::EdgeIndex, Self::OptionalEdgeIndex>;

    /// Returns true if this graph contains the given node index.
    fn contains_node_index(&self, node_id: Self::NodeIndex) -> bool;

    /// Returns true if this graph contains the given edge index.
    fn contains_edge_index(&self, edge_id: Self::EdgeIndex) -> bool;

    /// Returns the amount of nodes in this graph.
    fn node_count(&self) -> usize;

    /// Returns the amount of edges in this graph.
    fn edge_count(&self) -> usize;

    /// Returns a reference to the node data associated with the given node id, or None if there is no such node.
    fn node_data(&self, node_id: Self::NodeIndex) -> &Self::NodeData;

    /// Returns a reference to the edge data associated with the given edge id, or None if there is no such edge.
    fn edge_data(&self, edge_id: Self::EdgeIndex) -> &Self::EdgeData;

    /// Returns a mutable reference to the node data associated with the given node id, or None if there is no such node.
    fn node_data_mut(&mut self, node_id: Self::NodeIndex) -> &mut Self::NodeData;

    /// Returns a mutable reference to the edge data associated with the given edge id, or None if there is no such edge.
    fn edge_data_mut(&mut self, edge_id: Self::EdgeIndex) -> &mut Self::EdgeData;

    /// Returns true if the graph contains an edge `(from, to)`.
    fn contains_edge_between(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> bool;

    /// Returns the amount of edges `(from, to)`.
    fn edge_count_between(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> usize;

    /// Returns the endpoints of an edge.
    fn edge_endpoints(&self, edge_id: Self::EdgeIndex) -> Edge<Self::NodeIndex>;

    /// Returns true if the graph is empty, i.e. contains no nodes or edges.
    fn is_empty(&self) -> bool {
        // Zero nodes must imply zero edges.
        debug_assert!(self.node_count() != 0 || self.edge_count() == 0);
        self.node_count() == 0
    }
}

/// A container that allows adding and removing nodes and edges.
pub trait MutableGraphContainer: ImmutableGraphContainer {
    /// Adds a new node with the given `NodeData` to the graph.
    fn add_node(&mut self, node_data: Self::NodeData) -> Self::NodeIndex;

    /// Adds a new edge with the given `EdgeData` to the graph.
    fn add_edge(
        &mut self,
        from: Self::NodeIndex,
        to: Self::NodeIndex,
        edge_data: Self::EdgeData,
    ) -> Self::EdgeIndex;

    /// Removes the node with the given id from the graph.
    /// Note that this may change the ids of existing nodes.
    fn remove_node(&mut self, node_id: Self::NodeIndex) -> Option<Self::NodeData>;

    /// Removes the edge with the given id from the graph.
    /// Note that this may change the ids of existing edges.
    fn remove_edge(&mut self, edge_id: Self::EdgeIndex) -> Option<Self::EdgeData>;

    /// Removes the edges with the given ids from the graph.
    /// The ids are expected to be given in sorted order.
    ///
    /// Note that this may change the ids of existing edges.
    fn remove_edges_sorted(&mut self, edge_ids: &[Self::EdgeIndex]);

    /// Removes all nodes and edges from the graph.
    fn clear(&mut self);
}

/// A graph that can be navigated, i.e. that can iterate the neighbors of its nodes.
pub trait NavigableGraph<'a>: ImmutableGraphContainer + Sized {
    /// The iterator type used to iterate over the outgoing neighbors of a node.
    type OutNeighbors: Iterator<Item = Neighbor<Self::NodeIndex, Self::EdgeIndex>>;
    /// The iterator type used to iterate over the incoming neighbors of a node.
    type InNeighbors: Iterator<Item = Neighbor<Self::NodeIndex, Self::EdgeIndex>>;
    /// The iterator type used to iterate over the edges between to nodes.
    type EdgesBetween: Iterator<Item = Self::EdgeIndex>;

    /// Returns an iterator over the outgoing neighbors of the given node.
    fn out_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::OutNeighbors;
    /// Returns an iterator over the incoming neighbors of the given node.
    fn in_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::InNeighbors;

    /// Returns an iterator over the edges `(from_node_id, to_node_id)`.
    fn edges_between(
        &'a self,
        from_node_id: Self::NodeIndex,
        to_node_id: Self::NodeIndex,
    ) -> Self::EdgesBetween;

    /// Returns the amount of outgoing edges from a node.
    fn out_degree(&'a self, node_id: Self::NodeIndex) -> usize {
        self.out_neighbors(node_id).count()
    }

    /// Returns the amount of incoming edges to a node.
    fn in_degree(&'a self, node_id: Self::NodeIndex) -> usize {
        self.in_neighbors(node_id).count()
    }

    /// Returns true if the given node has indegree == 1 and outdegree == 1.
    fn is_biunivocal_node(&'a self, node_id: Self::NodeIndex) -> bool {
        self.in_degree(node_id) == 1 && self.out_degree(node_id) == 1
    }

    /// Returns true if the given node has indegree > 1 and outdegree > 1.
    fn is_bivalent_node(&'a self, node_id: Self::NodeIndex) -> bool {
        self.in_degree(node_id) > 1 && self.out_degree(node_id) > 1
    }

    /// Returns true if the given edge's tail has outdegree > 1.
    fn is_split_edge(&'a self, edge_id: Self::EdgeIndex) -> bool {
        self.out_degree(self.edge_endpoints(edge_id).from_node) > 1
    }

    /// Returns true if the given edge's head has indegree > 1.
    fn is_join_edge(&'a self, edge_id: Self::EdgeIndex) -> bool {
        self.in_degree(self.edge_endpoints(edge_id).to_node) > 1
    }

    /// Returns true if the given node has outdegree > 1.
    fn is_split_node(&'a self, node_id: Self::NodeIndex) -> bool {
        self.out_degree(node_id) > 1
    }

    /// Returns true if the given node has indegree > 1.
    fn is_join_node(&'a self, node_id: Self::NodeIndex) -> bool {
        self.in_degree(node_id) > 1
    }
}

/// A helper trait to get the correct walk type from a graph.
/// This is the factory pattern, where a graph is a factory for walks.
pub trait WalkableGraph: GraphBase + Sized {
    /// Create a node-centric walk over the given nodes in this graph.
    fn create_node_walk<
        WalkType: for<'a> NodeWalk<'a, Self, SubwalkType> + FromIterator<Self::NodeIndex>,
        SubwalkType: for<'a> NodeWalk<'a, Self, SubwalkType> + ?Sized,
    >(
        &self,
        walk: &[Self::NodeIndex],
    ) -> WalkType {
        walk.iter().copied().collect()
    }

    /// Create an empty node-centric walk in this graph.
    fn create_empty_node_walk<
        WalkType: for<'a> NodeWalk<'a, Self, SubwalkType> + Default,
        SubwalkType: for<'a> NodeWalk<'a, Self, SubwalkType> + ?Sized,
    >(
        &self,
    ) -> WalkType {
        Default::default()
    }

    /// Create an edge-centric walk over the given edges in this graph.
    fn create_edge_walk<
        WalkType: for<'a> EdgeWalk<'a, Self, SubwalkType> + FromIterator<Self::EdgeIndex>,
        SubwalkType: for<'a> EdgeWalk<'a, Self, SubwalkType> + ?Sized,
    >(
        &self,
        walk: &[Self::EdgeIndex],
    ) -> WalkType {
        walk.iter().copied().collect()
    }

    /// Create an empty edge-centric walk in this graph.
    fn create_empty_edge_walk<
        WalkType: for<'a> EdgeWalk<'a, Self, SubwalkType> + Default,
        SubwalkType: for<'a> EdgeWalk<'a, Self, SubwalkType> + ?Sized,
    >(
        &self,
    ) -> WalkType {
        Default::default()
    }
}
impl<Graph: GraphBase> WalkableGraph for Graph {}

/// A graph implementing all common graph traits that do not require mutable access.
/// This is a useful shortcut for generic type bounds when the graph should not be mutated.
pub trait StaticGraph:
    ImmutableGraphContainer + for<'a> NavigableGraph<'a> + WalkableGraph
{
}
impl<T: ImmutableGraphContainer + for<'a> NavigableGraph<'a> + WalkableGraph> StaticGraph for T {}

/// A graph implementing all common graph traits, including those requiring mutable access.
/// This is a useful shortcut for generic type bounds when the graph should be mutated.
pub trait DynamicGraph: StaticGraph + MutableGraphContainer {}
impl<T: StaticGraph + MutableGraphContainer> DynamicGraph for T {}

/// An edge represented as a pair of node indices.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Edge<NodeIndex> {
    /// The tail of this edge.
    pub from_node: NodeIndex,
    /// The head of this edge.
    pub to_node: NodeIndex,
}

/// The neighbor of a node, given as the edge used to reach the neighbor node as well as the neighbor node itself.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Neighbor<NodeIndex, EdgeIndex> {
    /// The edge used to reach the neighboring node.
    pub edge_id: EdgeIndex,
    /// The neighboring node.
    pub node_id: NodeIndex,
}

/// An enum encoding an index that can either be a node index or an edge index.
#[derive(Debug, Eq, PartialEq, Clone)]
pub enum NodeOrEdge<NodeIndex, EdgeIndex> {
    /// A node index.
    Node(NodeIndex),
    /// An edge index.
    Edge(EdgeIndex),
}
