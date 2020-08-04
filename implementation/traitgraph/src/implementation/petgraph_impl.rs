use crate::index::{GraphIndex, GraphIndices};
use crate::interface::{
    DynamicGraph, GraphBase, ImmutableGraphContainer, MutableGraphContainer, NavigableGraph,
    Neighbor,
};
use num_traits::{PrimInt, ToPrimitive};
use petgraph::graph::{DiGraph, Edges};
use petgraph::visit::EdgeRef;
use petgraph::{Directed, Direction};
use std::iter::Map;

pub use petgraph;

pub fn new<NodeData: 'static + Clone, EdgeData: 'static + Clone>(
) -> impl DynamicGraph<NodeData = NodeData, EdgeData = EdgeData> + Default + Clone {
    DiGraph::<NodeData, EdgeData, usize>::default()
}

impl<NodeData, EdgeData> GraphBase for DiGraph<NodeData, EdgeData, usize> {
    type NodeData = NodeData;
    type EdgeData = EdgeData;
    type OptionalNodeIndex = crate::index::OptionalNodeIndex<usize>;
    type OptionalEdgeIndex = crate::index::OptionalEdgeIndex<usize>;
    type NodeIndex = crate::index::NodeIndex<usize>;
    type EdgeIndex = crate::index::EdgeIndex<usize>;
}

impl<NodeData, EdgeData> ImmutableGraphContainer for DiGraph<NodeData, EdgeData, usize> {
    fn node_indices(&self) -> GraphIndices<Self::NodeIndex, Self::OptionalNodeIndex> {
        GraphIndices::from((0, self.node_count()))
    }

    fn edge_indices(&self) -> GraphIndices<Self::EdgeIndex, Self::OptionalEdgeIndex> {
        GraphIndices::from((0, self.edge_count()))
    }

    fn contains_node_index(&self, node_id: Self::NodeIndex) -> bool {
        self.node_weight(node_id.into()).is_some()
    }

    fn contains_edge_index(&self, edge_id: Self::EdgeIndex) -> bool {
        self.edge_weight(edge_id.into()).is_some()
    }

    fn node_count(&self) -> usize {
        self.node_count()
    }

    fn edge_count(&self) -> usize {
        self.edge_count()
    }

    fn node_data(&self, node_id: Self::NodeIndex) -> &Self::NodeData {
        self.node_weight(node_id.into()).unwrap()
    }

    fn edge_data(&self, edge_id: Self::EdgeIndex) -> &Self::EdgeData {
        self.edge_weight(edge_id.into()).unwrap()
    }

    fn node_data_mut(&mut self, node_id: Self::NodeIndex) -> &mut Self::NodeData {
        self.node_weight_mut(node_id.into()).unwrap()
    }

    fn edge_data_mut(&mut self, edge_id: Self::EdgeIndex) -> &mut Self::EdgeData {
        self.edge_weight_mut(edge_id.into()).unwrap()
    }

    fn contains_edge(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> bool {
        self.edges_connecting(from.into(), to.into())
            .next()
            .is_some()
    }
}

impl<NodeData, EdgeData> MutableGraphContainer for DiGraph<NodeData, EdgeData, usize> {
    fn add_node(&mut self, node_data: NodeData) -> Self::NodeIndex {
        self.add_node(node_data).index().into()
    }

    fn add_edge(
        &mut self,
        from: Self::NodeIndex,
        to: Self::NodeIndex,
        edge_data: EdgeData,
    ) -> Self::EdgeIndex {
        self.add_edge(from.into(), to.into(), edge_data)
            .index()
            .into()
    }

    fn remove_node(&mut self, node_id: Self::NodeIndex) -> Option<NodeData> {
        self.remove_node(node_id.into())
    }

    fn remove_edge(&mut self, edge_id: Self::EdgeIndex) -> Option<EdgeData> {
        self.remove_edge(edge_id.into())
    }
}

type PetgraphNeighborTranslator<'a, EdgeData, NodeIndex, EdgeIndex> = Map<
    Edges<'a, EdgeData, Directed, usize>,
    fn(petgraph::graph::EdgeReference<'a, EdgeData, usize>) -> Neighbor<NodeIndex, EdgeIndex>,
>;

impl<'a, NodeData, EdgeData: 'a> NavigableGraph<'a> for DiGraph<NodeData, EdgeData, usize> {
    type OutNeighbors = PetgraphNeighborTranslator<
        'a,
        EdgeData,
        <Self as GraphBase>::NodeIndex,
        <Self as GraphBase>::EdgeIndex,
    >;
    type InNeighbors = PetgraphNeighborTranslator<
        'a,
        EdgeData,
        <Self as GraphBase>::NodeIndex,
        <Self as GraphBase>::EdgeIndex,
    >;

    fn out_neighbors(&'a self, node_id: <Self as GraphBase>::NodeIndex) -> Self::OutNeighbors where
    {
        debug_assert!(node_id < self.node_count().into());
        self.edges_directed(node_id.into(), Direction::Outgoing)
            .map(|edge| Neighbor {
                edge_id: <Self as GraphBase>::EdgeIndex::from(edge.id().index()),
                node_id: <Self as GraphBase>::NodeIndex::from(edge.target().index()),
            })
    }

    fn in_neighbors(&'a self, node_id: <Self as GraphBase>::NodeIndex) -> Self::InNeighbors {
        debug_assert!(node_id < self.node_count().into());
        self.edges_directed(node_id.into(), Direction::Incoming)
            .map(|edge| Neighbor {
                edge_id: <Self as GraphBase>::EdgeIndex::from(edge.id().index()),
                node_id: <Self as GraphBase>::NodeIndex::from(edge.source().index()),
            })
    }
}

impl<IndexType: PrimInt + ToPrimitive + petgraph::graph::IndexType>
    From<crate::index::NodeIndex<IndexType>> for petgraph::graph::NodeIndex<IndexType>
{
    fn from(index: crate::index::NodeIndex<IndexType>) -> Self {
        petgraph::graph::NodeIndex::new(index.as_usize())
    }
}

impl<IndexType: PrimInt + ToPrimitive + petgraph::graph::IndexType>
    From<crate::index::EdgeIndex<IndexType>> for petgraph::graph::EdgeIndex<IndexType>
{
    fn from(index: crate::index::EdgeIndex<IndexType>) -> Self {
        petgraph::graph::EdgeIndex::new(index.as_usize())
    }
}
