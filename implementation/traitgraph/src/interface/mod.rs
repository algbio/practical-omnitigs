use crate::index::{GraphIndex, GraphIndices, OptionalGraphIndex};
use crate::walks::{EdgeWalk, NodeWalk};

/// Contains the associated types of a graph.
pub trait GraphBase {
    type NodeData;
    type EdgeData;
    type OptionalNodeIndex: OptionalGraphIndex<Self::NodeIndex>;
    type OptionalEdgeIndex: OptionalGraphIndex<Self::EdgeIndex>;
    type NodeIndex: GraphIndex<Self::OptionalNodeIndex>;
    type EdgeIndex: GraphIndex<Self::OptionalEdgeIndex>;
}

pub trait ImmutableGraphContainer: GraphBase {
    fn node_indices(&self) -> GraphIndices<Self::NodeIndex, Self::OptionalNodeIndex>;

    fn edge_indices(&self) -> GraphIndices<Self::EdgeIndex, Self::OptionalEdgeIndex>;

    fn contains_node_index(&self, node_id: Self::NodeIndex) -> bool;

    fn contains_edge_index(&self, edge_id: Self::EdgeIndex) -> bool;

    fn node_count(&self) -> usize;

    fn edge_count(&self) -> usize;

    /// Returns the node data associated with the given node id, or None if there is no such node.
    fn node_data(&self, node_id: Self::NodeIndex) -> &Self::NodeData;

    /// Returns the edge data associated with the given edge id, or None if there is no such edge.
    fn edge_data(&self, edge_id: Self::EdgeIndex) -> &Self::EdgeData;

    fn node_data_mut(&mut self, node_id: Self::NodeIndex) -> &mut Self::NodeData;

    fn edge_data_mut(&mut self, edge_id: Self::EdgeIndex) -> &mut Self::EdgeData;

    fn contains_edge(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> bool;

    fn edge_count_between(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> usize;

    fn edge_endpoints(&self, edge_id: Self::EdgeIndex) -> Edge<Self::NodeIndex>;

    /// Returns true if the graph is empty, i.e. contains no nodes or edges.
    fn is_empty(&self) -> bool {
        // Zero nodes must imply zero edges.
        debug_assert!(self.node_count() != 0 || self.edge_count() == 0);
        self.node_count() == 0
    }
}

pub trait MutableGraphContainer: GraphBase {
    fn add_node(&mut self, node_data: Self::NodeData) -> Self::NodeIndex;

    fn add_edge(
        &mut self,
        from: Self::NodeIndex,
        to: Self::NodeIndex,
        edge_data: Self::EdgeData,
    ) -> Self::EdgeIndex;

    fn remove_node(&mut self, node_id: Self::NodeIndex) -> Option<Self::NodeData>;

    fn remove_edge(&mut self, edge_id: Self::EdgeIndex) -> Option<Self::EdgeData>;
}

pub trait NavigableGraph<'a>: GraphBase {
    type OutNeighbors: Iterator<Item = Neighbor<Self::NodeIndex, Self::EdgeIndex>>;
    type InNeighbors: Iterator<Item = Neighbor<Self::NodeIndex, Self::EdgeIndex>>;
    type OutNeighborsTo: Iterator<Item = Self::EdgeIndex>;
    type InNeighborsFrom: Iterator<Item = Self::EdgeIndex>;

    fn out_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::OutNeighbors;
    fn in_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::InNeighbors;

    fn out_neighbors_to(
        &'a self,
        from_node_id: Self::NodeIndex,
        to_node_id: Self::NodeIndex,
    ) -> Self::OutNeighborsTo;
    fn in_neighbors_from(
        &'a self,
        to_node_id: Self::NodeIndex,
        from_node_id: Self::NodeIndex,
    ) -> Self::InNeighborsFrom;

    fn out_degree(&'a self, node_id: Self::NodeIndex) -> usize {
        self.out_neighbors(node_id).count()
    }

    fn in_degree(&'a self, node_id: Self::NodeIndex) -> usize {
        self.in_neighbors(node_id).count()
    }

    /// Returns true if the given node has indegree == outdegree == 1.
    fn is_biunivocal_node(&'a self, node_id: Self::NodeIndex) -> bool {
        self.in_degree(node_id) == 1 && self.out_degree(node_id) == 1
    }

    /// Returns true if the given node has indegree > 1 and outdegree > 1.
    fn is_bivalent_node(&'a self, node_id: Self::NodeIndex) -> bool {
        self.in_degree(node_id) > 1 && self.out_degree(node_id) > 1
    }
}

/// A helper trait to get the correct walk type from a graph.
/// This is the factory pattern, where a graph is a factory for walks.
pub trait WalkableGraph: GraphBase + Sized {
    fn create_node_walk<WalkType: for<'a> NodeWalk<'a, Self>>(
        &self,
        walk: &[Self::NodeIndex],
    ) -> WalkType {
        WalkType::from(&walk)
    }

    fn create_empty_node_walk<WalkType: for<'a> NodeWalk<'a, Self>>(&self) -> WalkType {
        self.create_node_walk(&[])
    }
    fn create_edge_walk<WalkType: for<'a> EdgeWalk<'a, Self>>(
        &self,
        walk: &[Self::EdgeIndex],
    ) -> WalkType {
        WalkType::from(&walk)
    }

    fn create_empty_edge_walk<WalkType: for<'a> EdgeWalk<'a, Self>>(&self) -> WalkType {
        self.create_edge_walk(&[])
    }
}
impl<Graph: GraphBase> WalkableGraph for Graph {}

pub trait StaticGraph:
    ImmutableGraphContainer + for<'a> NavigableGraph<'a> + WalkableGraph
{
}
impl<T: ImmutableGraphContainer + for<'a> NavigableGraph<'a> + WalkableGraph> StaticGraph for T {}

pub trait DynamicGraph: StaticGraph + MutableGraphContainer {}
impl<T: StaticGraph + MutableGraphContainer> DynamicGraph for T {}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Edge<NodeIndex> {
    pub from_node: NodeIndex,
    pub to_node: NodeIndex,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Neighbor<NodeIndex, EdgeIndex> {
    pub edge_id: EdgeIndex,
    pub node_id: NodeIndex,
}
