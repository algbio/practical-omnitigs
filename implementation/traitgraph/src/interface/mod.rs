use crate::{EdgeIndex, EdgeIndices, NodeIndex, NodeIndices};
use num_traits::PrimInt;

pub trait ImmutableGraphContainer {
    type NodeData;
    type EdgeData;
    type IndexType;

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

pub trait MutableGraphContainer<NodeData, EdgeData, IndexType> {
    fn add_node(&mut self, node_data: NodeData) -> NodeIndex<IndexType>;

    fn add_edge(
        &mut self,
        from: NodeIndex<IndexType>,
        to: NodeIndex<IndexType>,
        edge_data: EdgeData,
    ) -> EdgeIndex<IndexType>;

    fn remove_node(&mut self, node_id: NodeIndex<IndexType>) -> Option<NodeData>;

    fn remove_edge(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<EdgeData>;
}

pub trait NavigableGraph<'a, NodeData, EdgeData, IndexType> {
    type OutNeighbors: IntoIterator<Item = Neighbor<IndexType>>;
    type InNeighbors: IntoIterator<Item = Neighbor<IndexType>>;

    fn out_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::OutNeighbors>;

    fn in_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::InNeighbors>;
}

pub trait StaticGraph<NodeData, EdgeData, IndexType: PrimInt>:
    ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
    + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType>
{
}
impl<
        NodeData,
        EdgeData,
        IndexType: PrimInt,
        T: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
            + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType>,
    > StaticGraph<NodeData, EdgeData, IndexType> for T
{
}

pub trait DynamicGraph<NodeData, EdgeData, IndexType: PrimInt>:
    ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
    + MutableGraphContainer<NodeData, EdgeData, IndexType>
    + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType>
{
}
impl<
        NodeData,
        EdgeData,
        IndexType: PrimInt,
        T: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>
            + MutableGraphContainer<NodeData, EdgeData, IndexType>
            + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType>,
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
