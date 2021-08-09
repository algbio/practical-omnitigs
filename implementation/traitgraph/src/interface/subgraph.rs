use crate::interface::{
    Edge, GraphBase, GraphIndices, ImmutableGraphContainer, NavigableGraph, Neighbor,
};
use std::borrow::Borrow;

/// A subset of nodes and edges of a graph.
/// There is not restriction on which nodes or edges must be contained in combination, so e.g. an edge from node 1 to node 4 can be contained, even if node 1 and node 4 are both not contained.
pub trait DecoratingSubgraph {
    /// The type of the associated parent graph.
    type ParentGraph: GraphBase;
    /// The type of the reference to the associated parent graph.
    type ParentGraphRef: Borrow<Self::ParentGraph>;

    /// Constructs a subgraph from the given graph without any nodes or edges.
    /// If not defined otherwise in the implementation, all node and edge ids of the given graph are valid arguments for the methods of this trait on the returned object.
    fn new_empty(graph: Self::ParentGraphRef) -> Self;

    /// Constructs a subgraph from the given graph with all nodes and edges.
    /// If not defined otherwise in the implementation, all node and edge ids of the given graph are valid arguments for the methods of this trait on the returned object.
    fn new_full(graph: Self::ParentGraphRef) -> Self;

    /// Removes all nodes and edges from the subgraph.
    fn clear(&mut self);

    /// Returns a reference to the original graph.
    fn parent_graph(&self) -> &Self::ParentGraph;

    /// Returns true if the given node id is part of the subgraph.
    fn contains_node(&self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex) -> bool;

    /// Returns true if the given edge id is part of the subgraph.
    fn contains_edge(&self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex) -> bool;

    /// Adds the given node id to the subgraph.
    /// Note that some implementations may require an initial set of nodes to be known when a subgraph is created, and inserting nodes outside of this initial set might panic.
    fn add_node(&mut self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex);

    /// Adds the given edge id to the subgraph.
    /// Note that some implementations may require an initial set of nodes to be known when a subgraph is created, and inserting nodes outside of this initial set might panic.
    fn add_edge(&mut self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex);

    /// Removes the given node id from the graph.
    fn remove_node(&mut self, node_index: <Self::ParentGraph as GraphBase>::NodeIndex);

    /// Removes the given edge id from the graph.
    fn remove_edge(&mut self, edge_index: <Self::ParentGraph as GraphBase>::EdgeIndex);

    /// Returns the amount of nodes in the subgraph.
    fn node_count(&self) -> usize;

    /// Returns the amount of edges in the subgraph.
    fn edge_count(&self) -> usize;
}

impl<T: DecoratingSubgraph> GraphBase for T {
    type NodeData = <<Self as DecoratingSubgraph>::ParentGraph as GraphBase>::NodeData;
    type EdgeData = <<Self as DecoratingSubgraph>::ParentGraph as GraphBase>::EdgeData;
    type OptionalNodeIndex =
        <<Self as DecoratingSubgraph>::ParentGraph as GraphBase>::OptionalNodeIndex;
    type OptionalEdgeIndex =
        <<Self as DecoratingSubgraph>::ParentGraph as GraphBase>::OptionalEdgeIndex;
    type NodeIndex = <<Self as DecoratingSubgraph>::ParentGraph as GraphBase>::NodeIndex;
    type EdgeIndex = <<Self as DecoratingSubgraph>::ParentGraph as GraphBase>::EdgeIndex;
}

impl<T: DecoratingSubgraph> ImmutableGraphContainer for T
where
    T::ParentGraph: ImmutableGraphContainer + for<'a> NavigableGraph<'a>,
{
    fn node_indices(&self) -> GraphIndices<Self::NodeIndex, Self::OptionalNodeIndex> {
        unimplemented!("Will not implement if not necessary");
    }

    fn edge_indices(&self) -> GraphIndices<Self::EdgeIndex, Self::OptionalEdgeIndex> {
        unimplemented!("Will not implement if not necessary");
    }

    fn contains_node_index(&self, node_id: Self::NodeIndex) -> bool {
        <Self as DecoratingSubgraph>::contains_node(self, node_id)
    }

    fn contains_edge_index(&self, edge_id: Self::EdgeIndex) -> bool {
        <Self as DecoratingSubgraph>::contains_edge(self, edge_id)
    }

    fn node_count(&self) -> usize {
        <Self as DecoratingSubgraph>::node_count(self)
    }

    fn edge_count(&self) -> usize {
        <Self as DecoratingSubgraph>::edge_count(self)
    }

    fn node_data(&self, node_id: Self::NodeIndex) -> &Self::NodeData {
        self.parent_graph().node_data(node_id)
    }

    fn edge_data(&self, edge_id: Self::EdgeIndex) -> &Self::EdgeData {
        self.parent_graph().edge_data(edge_id)
    }

    fn node_data_mut(&mut self, _node_id: Self::NodeIndex) -> &mut Self::NodeData {
        unimplemented!("Cannot access parent graph mutably")
    }

    fn edge_data_mut(&mut self, _edge_id: Self::EdgeIndex) -> &mut Self::EdgeData {
        unimplemented!("Cannot access parent graph mutably")
    }

    fn contains_edge_between(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> bool {
        self.edge_count_between(from, to) > 0
    }

    fn edge_count_between(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> usize {
        self.parent_graph()
            .edges_between(from, to)
            .filter(|e| <Self as DecoratingSubgraph>::contains_edge(self, *e))
            .count()
    }

    fn edge_endpoints(&self, edge_id: Self::EdgeIndex) -> Edge<Self::NodeIndex> {
        self.parent_graph().edge_endpoints(edge_id)
    }
}

impl<'a, T: 'a + DecoratingSubgraph> NavigableGraph<'a> for T
where
    T::ParentGraph: ImmutableGraphContainer + for<'b> NavigableGraph<'b>,
{
    //type OutNeighbors = <<Self as DecoratingSubgraph>::ParentGraph as NavigableGraph<'a>>::OutNeighbors;//std::iter::Filter<<<Self as DecoratingSubgraph>::ParentGraph as NavigableGraph<'a>>::OutNeighbors, fn(&Neighbor<<Self as GraphBase>::NodeIndex,<Self as GraphBase>::EdgeIndex>)->bool>;
    type OutNeighbors = EdgeFilteredNeighborIterator<
        'a,
        T,
        <<Self as DecoratingSubgraph>::ParentGraph as NavigableGraph<'a>>::OutNeighbors,
    >;
    type InNeighbors = EdgeFilteredNeighborIterator<
        'a,
        T,
        <<Self as DecoratingSubgraph>::ParentGraph as NavigableGraph<'a>>::InNeighbors,
    >;
    type EdgesBetween = EdgeFilteredEdgeIterator<
        'a,
        T,
        <<Self as DecoratingSubgraph>::ParentGraph as NavigableGraph<'a>>::EdgesBetween,
    >;

    fn out_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::OutNeighbors {
        EdgeFilteredNeighborIterator {
            graph: self,
            iter: self.parent_graph().out_neighbors(node_id),
        }
    }

    fn in_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::InNeighbors {
        EdgeFilteredNeighborIterator {
            graph: self,
            iter: self.parent_graph().in_neighbors(node_id),
        }
    }
    fn edges_between(
        &'a self,
        from_node_id: Self::NodeIndex,
        to_node_id: Self::NodeIndex,
    ) -> Self::EdgesBetween {
        EdgeFilteredEdgeIterator {
            graph: self,
            iter: self.parent_graph().edges_between(from_node_id, to_node_id),
        }
    }
}

/// An iterator adapter that filters neighbors that are not part of a subgraph.
/// A neighbor is considered to not be part of a subgraph if its connecting edge is not part of the subgraph.
/// The neighboring node does not influence filtering.
pub struct EdgeFilteredNeighborIterator<'a, Graph, SourceIterator> {
    graph: &'a Graph,
    iter: SourceIterator,
}

impl<
        'a,
        Graph: DecoratingSubgraph,
        SourceIterator: Iterator<Item = Neighbor<<Graph as GraphBase>::NodeIndex, <Graph as GraphBase>::EdgeIndex>>,
    > Iterator for EdgeFilteredNeighborIterator<'a, Graph, SourceIterator>
{
    type Item = Neighbor<<Graph as GraphBase>::NodeIndex, <Graph as GraphBase>::EdgeIndex>;

    fn next(&mut self) -> Option<Self::Item> {
        for neighbor in &mut self.iter {
            if self.graph.contains_edge(neighbor.edge_id) {
                return Some(neighbor);
            }
        }

        None
    }
}

/// An iterator adapter that filters edges that are not part of a subgraph.
pub struct EdgeFilteredEdgeIterator<'a, Graph, SourceIterator> {
    graph: &'a Graph,
    iter: SourceIterator,
}

impl<
        'a,
        Graph: DecoratingSubgraph,
        SourceIterator: Iterator<Item = <Graph as GraphBase>::EdgeIndex>,
    > Iterator for EdgeFilteredEdgeIterator<'a, Graph, SourceIterator>
{
    type Item = <Graph as GraphBase>::EdgeIndex;

    fn next(&mut self) -> Option<Self::Item> {
        for edge in &mut self.iter {
            if self.graph.contains_edge(edge) {
                return Some(edge);
            }
        }

        None
    }
}
