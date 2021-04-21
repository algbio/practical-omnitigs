use crate::interface::dynamic_bigraph::DynamicEdgeCentricBigraph;
use crate::interface::dynamic_bigraph::DynamicNodeCentricBigraph;
use crate::interface::{
    dynamic_bigraph::DynamicBigraph,
    static_bigraph::StaticBigraph,
    static_bigraph::{
        StaticBigraphFromDigraph, StaticEdgeCentricBigraph, StaticNodeCentricBigraph,
    },
    BidirectedData,
};
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use traitgraph::index::{GraphIndex, GraphIndices, OptionalGraphIndex};
use traitgraph::interface::{
    DynamicGraph, Edge, GraphBase, ImmutableGraphContainer, MutableGraphContainer, NavigableGraph,
    StaticGraph,
};

/// Represent arbitrary bigraphs with petgraph.
pub type PetBigraph<NodeData, EdgeData> =
    crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            NodeData,
            EdgeData,
            usize,
        >,
    >;

/// Wrapper for a static graph that adds a mirror node mapping function.
///
/// Bigraphs can be represented with this struct by creating their topology as normal directed graph where each binode is split into its two parts.
/// The binode mapping function then associates the parts with each other.
///
/// ```rust
/// use bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
/// use bigraph::interface::static_bigraph::{StaticBigraph, StaticBigraphFromDigraph};
/// use bigraph::traitgraph::interface::MutableGraphContainer;
/// use bigraph::traitgraph::implementation::petgraph_impl;
/// use bigraph::interface::BidirectedData;
///
/// #[derive(Clone, Eq, PartialEq, Hash, Debug)]
/// struct NodeData(i32);
/// impl BidirectedData for NodeData {
///     fn mirror(&self) -> Self {
///         Self(1000 - self.0)
///     }
/// }
///
/// let mut graph = petgraph_impl::new();
/// let n1 = graph.add_node(NodeData(0));
/// let n2 = graph.add_node(NodeData(1000));
/// graph.add_edge(n1.clone(), n2.clone(), ());
/// graph.add_edge(n2.clone(), n1.clone(), ());
/// let bigraph = NodeBigraphWrapper::new(graph);
/// assert_eq!(Some(n2.clone()), bigraph.mirror_node(n1.clone()));
/// assert_eq!(Some(n1.clone()), bigraph.mirror_node(n2.clone()));
/// ```
#[derive(Debug, Clone)]
pub struct NodeBigraphWrapper<Topology: GraphBase> {
    /// The underlying topology of the bigraph.
    pub topology: Topology,
    binode_map: Vec<Topology::OptionalNodeIndex>,
}

impl<'a, Topology: GraphBase> GraphBase for NodeBigraphWrapper<Topology> {
    type NodeData = Topology::NodeData;
    type EdgeData = Topology::EdgeData;
    type OptionalNodeIndex = Topology::OptionalNodeIndex;
    type OptionalEdgeIndex = Topology::OptionalEdgeIndex;
    type NodeIndex = Topology::NodeIndex;
    type EdgeIndex = Topology::EdgeIndex;
}

impl<Topology: StaticGraph> StaticBigraph for NodeBigraphWrapper<Topology> {
    fn mirror_node(&self, node_id: Self::NodeIndex) -> Option<Self::NodeIndex> {
        self.binode_map[node_id.as_usize()].into()
    }
}

impl<Topology: StaticGraph> NodeBigraphWrapper<Topology>
where
    <Self as GraphBase>::NodeData: BidirectedData + Eq + Hash + Debug,
{
    fn new_internal(topology: Topology, checked: bool) -> Self {
        //let mut data_map = HashMap::new();
        let mut data_map: HashMap<<Self as GraphBase>::NodeData, <Self as GraphBase>::NodeIndex> =
            HashMap::new();
        let mut binode_map =
            vec![<Self as GraphBase>::OptionalNodeIndex::new_none(); topology.node_count()];

        for node_index in topology.node_indices() {
            let node_data = topology.node_data(node_index);

            if let Some(mirror_index) = data_map.get(node_data).cloned() {
                //assert_ne!(node_index, mirror_index);
                assert!(!binode_map[node_index.as_usize()].is_valid());
                assert!(!binode_map[mirror_index.as_usize()].is_valid());
                assert_eq!(node_data, &topology.node_data(mirror_index).mirror());
                binode_map[node_index.as_usize()] = mirror_index.into();
                binode_map[mirror_index.as_usize()] = node_index.into();
                data_map.remove(node_data);
            } else {
                let mirror_data = node_data.mirror();
                //assert_ne!(&mirror_data, node_data);
                assert_eq!(None, data_map.insert(mirror_data, node_index));
            }
        }

        if checked {
            assert!(data_map.is_empty());
        } else {
            for node_index in topology.node_indices() {
                let node_data = topology.node_data(node_index);
                assert!(!data_map.contains_key(node_data));
            }
        }
        Self {
            topology,
            binode_map,
        }
    }
}

impl<Topology: StaticGraph> StaticBigraphFromDigraph for NodeBigraphWrapper<Topology>
where
    <Self as GraphBase>::NodeData: BidirectedData + Eq + Hash + Debug,
{
    type Topology = Topology;

    fn new(topology: Self::Topology) -> Self {
        Self::new_internal(topology, true)
    }

    fn new_unchecked(topology: Self::Topology) -> Self {
        Self::new_internal(topology, false)
    }
}

impl<Topology: DynamicGraph> DynamicBigraph for NodeBigraphWrapper<Topology> {
    fn set_mirror_nodes(&mut self, a: Self::NodeIndex, b: Self::NodeIndex) {
        assert!(self.contains_node_index(a));
        assert!(self.contains_node_index(b));
        self.binode_map[a.as_usize()] = b.into();
        self.binode_map[b.as_usize()] = a.into();
    }
}

impl<Topology: ImmutableGraphContainer> ImmutableGraphContainer for NodeBigraphWrapper<Topology> {
    fn node_indices(&self) -> GraphIndices<Self::NodeIndex, Self::OptionalNodeIndex> {
        self.topology.node_indices()
    }

    fn edge_indices(&self) -> GraphIndices<Self::EdgeIndex, Self::OptionalEdgeIndex> {
        self.topology.edge_indices()
    }

    fn contains_node_index(&self, node_index: Self::NodeIndex) -> bool {
        self.topology.contains_node_index(node_index)
    }

    fn contains_edge_index(&self, edge_index: Self::EdgeIndex) -> bool {
        self.topology.contains_edge_index(edge_index)
    }

    fn node_count(&self) -> usize {
        self.topology.node_count()
    }

    fn edge_count(&self) -> usize {
        self.topology.edge_count()
    }

    fn node_data(&self, node_id: Self::NodeIndex) -> &Self::NodeData {
        self.topology.node_data(node_id)
    }

    fn edge_data(&self, edge_id: Self::EdgeIndex) -> &Self::EdgeData {
        self.topology.edge_data(edge_id)
    }

    fn node_data_mut(&mut self, node_id: Self::NodeIndex) -> &mut Self::NodeData {
        self.topology.node_data_mut(node_id)
    }

    fn edge_data_mut(&mut self, edge_id: Self::EdgeIndex) -> &mut Self::EdgeData {
        self.topology.edge_data_mut(edge_id)
    }

    fn contains_edge_between(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> bool {
        self.topology.contains_edge_between(from, to)
    }

    fn edge_count_between(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> usize {
        self.topology.edge_count_between(from, to)
    }

    fn edge_endpoints(&self, edge_id: Self::EdgeIndex) -> Edge<Self::NodeIndex> {
        self.topology.edge_endpoints(edge_id)
    }
}

impl<Topology: MutableGraphContainer + StaticGraph> MutableGraphContainer
    for NodeBigraphWrapper<Topology>
{
    fn add_node(&mut self, node_data: Self::NodeData) -> Self::NodeIndex {
        self.binode_map.push(Self::OptionalNodeIndex::new_none());
        self.topology.add_node(node_data)
    }

    fn add_edge(
        &mut self,
        from: Self::NodeIndex,
        to: Self::NodeIndex,
        edge_data: Self::EdgeData,
    ) -> Self::EdgeIndex {
        self.topology.add_edge(from, to, edge_data)
    }

    fn remove_node(&mut self, node_id: Self::NodeIndex) -> Option<Self::NodeData> {
        self.topology.remove_node(node_id)
    }

    fn remove_edge(&mut self, edge_id: Self::EdgeIndex) -> Option<Self::EdgeData> {
        self.topology.remove_edge(edge_id)
    }

    fn remove_edges_sorted(&mut self, edge_ids: &[Self::EdgeIndex]) {
        self.topology.remove_edges_sorted(edge_ids)
    }

    fn clear(&mut self) {
        self.topology.clear();
        self.binode_map.clear();
    }
}

impl<'a, Topology: NavigableGraph<'a>> NavigableGraph<'a> for NodeBigraphWrapper<Topology> {
    type OutNeighbors = <Topology as NavigableGraph<'a>>::OutNeighbors;
    type InNeighbors = <Topology as NavigableGraph<'a>>::InNeighbors;
    type EdgesBetween = <Topology as NavigableGraph<'a>>::EdgesBetween;

    fn out_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::OutNeighbors {
        self.topology.out_neighbors(node_id)
    }

    fn in_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::InNeighbors {
        self.topology.in_neighbors(node_id)
    }

    fn edges_between(
        &'a self,
        from_node_id: Self::NodeIndex,
        to_node_id: Self::NodeIndex,
    ) -> Self::EdgesBetween {
        self.topology.edges_between(from_node_id, to_node_id)
    }
}

impl<Topology: Default + GraphBase> Default for NodeBigraphWrapper<Topology> {
    fn default() -> Self {
        Self {
            topology: Topology::default(),
            binode_map: vec![],
        }
    }
}

impl<Topology: StaticGraph> StaticNodeCentricBigraph for NodeBigraphWrapper<Topology> {}

impl<Topology: StaticGraph> StaticEdgeCentricBigraph for NodeBigraphWrapper<Topology> where
    <Topology as GraphBase>::EdgeData: BidirectedData + Eq
{
}

impl<Topology: DynamicGraph> DynamicNodeCentricBigraph for NodeBigraphWrapper<Topology>
where
    <Topology as GraphBase>::NodeData: BidirectedData,
    <Topology as GraphBase>::EdgeData: Clone,
{
}

impl<Topology: DynamicGraph> DynamicEdgeCentricBigraph for NodeBigraphWrapper<Topology> where
    <Topology as GraphBase>::EdgeData: BidirectedData + Eq
{
}

#[cfg(test)]
mod tests {
    use super::NodeBigraphWrapper;
    use crate::interface::dynamic_bigraph::DynamicNodeCentricBigraph;
    use crate::interface::{
        static_bigraph::StaticBigraph, static_bigraph::StaticBigraphFromDigraph, BidirectedData,
    };
    use crate::traitgraph::implementation::petgraph_impl;
    use crate::traitgraph::interface::{ImmutableGraphContainer, MutableGraphContainer};
    use traitgraph::index::OptionalGraphIndex;

    #[test]
    fn test_bigraph_creation() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        let n3 = graph.add_node(NodeData(2));
        let n4 = graph.add_node(NodeData(3));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        let bigraph = NodeBigraphWrapper::new(graph);

        assert_eq!(Some(n2), bigraph.mirror_node(n1));
        assert_eq!(Some(n1), bigraph.mirror_node(n2));
        assert_eq!(Some(n4), bigraph.mirror_node(n3));
        assert_eq!(Some(n3), bigraph.mirror_node(n4));
    }

    #[test]
    #[should_panic]
    fn test_bigraph_creation_unmapped_node() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(4));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new(graph);
    }

    #[test]
    #[should_panic]
    fn test_bigraph_creation_wrongly_mapped_node() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 == 4 {
                    3
                } else if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(4));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new(graph);
    }

    #[test]
    #[should_panic]
    fn test_bigraph_creation_self_mapped_node_without_mirror() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 == 4 {
                    4
                } else if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(4));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new(graph);
    }

    #[test]
    #[should_panic]
    fn test_bigraph_creation_self_mapped_node_with_mirror() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 == 4 {
                    4
                } else if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(4));
        graph.add_node(NodeData(5));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new(graph);
    }

    #[test]
    fn test_bigraph_unchecked_creation() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        let n3 = graph.add_node(NodeData(2));
        let n4 = graph.add_node(NodeData(3));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        let bigraph = NodeBigraphWrapper::new_unchecked(graph);

        assert_eq!(Some(n2), bigraph.mirror_node(n1));
        assert_eq!(Some(n1), bigraph.mirror_node(n2));
        assert_eq!(Some(n4), bigraph.mirror_node(n3));
        assert_eq!(Some(n3), bigraph.mirror_node(n4));
    }

    #[test]
    fn test_bigraph_unchecked_creation_unmapped_node() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        let n3 = graph.add_node(NodeData(2));
        let n4 = graph.add_node(NodeData(3));
        let n5 = graph.add_node(NodeData(4));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        let bigraph = NodeBigraphWrapper::new_unchecked(graph);

        assert_eq!(Some(n2), bigraph.mirror_node(n1));
        assert_eq!(Some(n1), bigraph.mirror_node(n2));
        assert_eq!(Some(n4), bigraph.mirror_node(n3));
        assert_eq!(Some(n3), bigraph.mirror_node(n4));
        assert_eq!(None, bigraph.mirror_node(n5));
    }

    #[test]
    #[should_panic]
    fn test_bigraph_unchecked_creation_wrongly_mapped_node_at_end() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 == 4 {
                    3
                } else if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(4));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new_unchecked(graph);
    }

    #[test]
    #[should_panic]
    fn test_bigraph_unchecked_creation_wrongly_mapped_node_at_beginning() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 == 4 {
                    3
                } else if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        graph.add_node(NodeData(4));
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new_unchecked(graph);
    }

    #[test]
    #[should_panic]
    fn test_bigraph_unchecked_creation_self_mapped_node_without_mirror() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 == 4 {
                    4
                } else if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(4));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new_unchecked(graph);
    }

    #[test]
    #[should_panic]
    fn test_bigraph_unchecked_creation_self_mapped_node_with_mirror() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 == 4 {
                    4
                } else if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(4));
        graph.add_node(NodeData(5));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new_unchecked(graph);
    }

    #[test]
    fn test_bigraph_verification() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        let bigraph = NodeBigraphWrapper::new(graph);
        assert!(bigraph.verify_node_pairing_without_self_mirrors());
        assert!(bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_verification_self_mapped_node() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        let mut bigraph = NodeBigraphWrapper::new(graph);
        bigraph.topology.add_node(NodeData(4));
        bigraph.binode_map.push(4usize.into());
        assert!(!bigraph.verify_node_pairing_without_self_mirrors());
        assert!(bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_verification_self_unmapped_node() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        let mut bigraph = NodeBigraphWrapper::new(graph);
        bigraph.topology.add_node(NodeData(4));
        bigraph.binode_map.push(OptionalGraphIndex::new_none());
        assert!(!bigraph.verify_node_pairing_without_self_mirrors());
        assert!(!bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_verification_wrongly_mapped_node() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(if self.0 % 2 == 0 {
                    self.0 + 1
                } else {
                    self.0 - 1
                })
            }
        }

        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_edge(n1, n2, ()); // Just to fix the EdgeData type parameter
        let mut bigraph = NodeBigraphWrapper::new(graph);
        bigraph.topology.add_node(NodeData(4));
        bigraph.binode_map.push(3usize.into());
        assert!(!bigraph.verify_node_pairing_without_self_mirrors());
        assert!(!bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_add_mirror_nodes() {
        #[derive(Eq, PartialEq, Debug, Hash, Clone)]
        struct NodeData(u32);
        impl BidirectedData for NodeData {
            fn mirror(&self) -> Self {
                Self(1000 - self.0)
            }
        }

        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(NodeData(0));
        let n1 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(997));
        graph.add_edge(n0, n1, ());
        let mut graph = NodeBigraphWrapper::new_unchecked(graph);
        assert!(!graph.verify_node_pairing());
        graph.add_mirror_nodes();
        assert!(graph.verify_node_pairing());
        assert_eq!(graph.node_count(), 8);
    }
}
