use crate::{EdgeIndex, EdgeIndices, NodeIndex, NodeIndices};

/// Contains the associated types of a graph.
pub trait GraphBase {
    type NodeData;
    type EdgeData;
    type IndexType;
}

pub trait ImmutableGraphContainer: GraphBase {
    fn node_indices(&self) -> NodeIndices<Self::IndexType>;

    fn edge_indices(&self) -> EdgeIndices<Self::IndexType>;

    fn contains_node_index(&self, node_id: NodeIndex<Self::IndexType>) -> bool {
        self.node_data(node_id).is_some()
    }

    fn contains_edge_index(&self, edge_id: EdgeIndex<Self::IndexType>) -> bool {
        self.edge_data(edge_id).is_some()
    }

    fn node_count(&self) -> usize;

    fn edge_count(&self) -> usize;

    /// Returns the node data associated with the given node id, or None if there is no such node.
    fn node_data(&self, node_id: NodeIndex<Self::IndexType>) -> Option<&Self::NodeData>;

    /// Returns the edge data associated with the given edge id, or None if there is no such edge.
    fn edge_data(&self, edge_id: EdgeIndex<Self::IndexType>) -> Option<&Self::EdgeData>;

    fn node_data_mut(&mut self, node_id: NodeIndex<Self::IndexType>)
        -> Option<&mut Self::NodeData>;

    fn edge_data_mut(&mut self, edge_id: EdgeIndex<Self::IndexType>)
        -> Option<&mut Self::EdgeData>;

    fn contains_edge(
        &self,
        from: NodeIndex<Self::IndexType>,
        to: NodeIndex<Self::IndexType>,
    ) -> bool;
}

pub trait MutableGraphContainer: GraphBase {
    fn add_node(&mut self, node_data: Self::NodeData) -> NodeIndex<Self::IndexType>;

    fn add_edge(
        &mut self,
        from: NodeIndex<Self::IndexType>,
        to: NodeIndex<Self::IndexType>,
        edge_data: Self::EdgeData,
    ) -> EdgeIndex<Self::IndexType>;

    fn remove_node(&mut self, node_id: NodeIndex<Self::IndexType>) -> Option<Self::NodeData>;

    fn remove_edge(&mut self, edge_id: EdgeIndex<Self::IndexType>) -> Option<Self::EdgeData>;
}

pub trait NavigableGraph<'a>: GraphBase {
    type OutNeighbors: IntoIterator<Item = Neighbor<Self::IndexType>>;
    type InNeighbors: IntoIterator<Item = Neighbor<Self::IndexType>>;

    fn out_neighbors(&'a self, node_id: NodeIndex<Self::IndexType>) -> Option<Self::OutNeighbors>;
    fn in_neighbors(&'a self, node_id: NodeIndex<Self::IndexType>) -> Option<Self::InNeighbors>;
}

pub trait StaticGraph: ImmutableGraphContainer + for<'a> NavigableGraph<'a> {}
impl<T: ImmutableGraphContainer + for<'a> NavigableGraph<'a>> StaticGraph for T {}

pub trait DynamicGraph<NodeData, EdgeData, IndexType>:
    ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
    + MutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
    + for<'a> NavigableGraph<'a, NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
{
}
impl<
        NodeData,
        EdgeData,
        IndexType,
        T: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
            + MutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
            + for<'a> NavigableGraph<
                'a,
                NodeData = NodeData,
                EdgeData = EdgeData,
                IndexType = IndexType,
            >,
    > DynamicGraph<NodeData, EdgeData, IndexType> for T
{
}

/*pub struct Edge<IndexType> {
    pub from: NodeIndex<IndexType>,
    pub to: NodeIndex<IndexType>,
}*/

pub struct Neighbor<IndexType> {
    pub edge_id: EdgeIndex<IndexType>,
    pub node_id: NodeIndex<IndexType>,
}
