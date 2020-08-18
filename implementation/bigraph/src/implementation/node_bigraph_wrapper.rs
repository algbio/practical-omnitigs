use crate::interface::{
    dynamic_bigraph::DynamicBigraph, static_bigraph::StaticBigraph,
    static_bigraph::StaticBigraphFromDigraph, BidirectedData,
};
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use traitgraph::index::{GraphIndex, GraphIndices, OptionalGraphIndex};
use traitgraph::interface::{
    DynamicGraph, Edge, GraphBase, ImmutableGraphContainer, MutableGraphContainer, NavigableGraph,
    StaticGraph,
};

/**
*   Wrapper for a static graph that adds a binode mapping function.
*
*   Bigraphs can be represented with this struct by creating their topology as normal directed graph where each binode is split into its two parts.
*   The binode mapping function then associates the parts with each other.
*
*   ```rust
*   use bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
*   use bigraph::interface::static_bigraph::{StaticBigraph, StaticBigraphFromDigraph};
*   use bigraph::traitgraph::interface::MutableGraphContainer;
*   use bigraph::traitgraph::implementation::petgraph_impl;
*
*   let mut graph = petgraph_impl::new();
*   let n1 = graph.add_node(0);
*   let n2 = graph.add_node(1);
*   graph.add_edge(n1.clone(), n2.clone(), ());
*   graph.add_edge(n2.clone(), n1.clone(), ());
*   let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 {n + 1} else {n - 1});
*   assert_eq!(Some(n2.clone()), bigraph.partner_node(n1.clone()));
*   assert_eq!(Some(n1.clone()), bigraph.partner_node(n2.clone()));
*   ```
*/
#[derive(Debug)]
pub struct NodeBigraphWrapper<Topology: GraphBase> {
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
    fn partner_node(&self, node_id: Self::NodeIndex) -> Option<Self::NodeIndex> {
        if let Some(partner_node) = self.binode_map.get(node_id.as_usize()).cloned() {
            if let Some(partner_node) = Into::<Option<Self::NodeIndex>>::into(partner_node) {
                Some(partner_node)
            } else {
                None
            }
        } else {
            None
        }
    }
}

impl<Topology: StaticGraph> NodeBigraphWrapper<Topology>
where
    <Self as GraphBase>::NodeData: Eq + Hash + Debug,
{
    fn new_internal(
        topology: Topology,
        binode_mapping_function: fn(
            &<Self as GraphBase>::NodeData,
        ) -> <Self as GraphBase>::NodeData,
        checked: bool,
    ) -> Self {
        let mut data_map = HashMap::new(); //: HashMap<NodeData, Self::NodeIndex> = HashMap::new();
        let mut binode_map =
            vec![<Self as GraphBase>::OptionalNodeIndex::new_none(); topology.node_count()];

        for node_index in topology.node_indices() {
            let node_data = topology.node_data(node_index);

            if let Some(partner_index) = data_map.get(node_data).cloned() {
                assert_ne!(node_index, partner_index);
                assert!(!binode_map[node_index.as_usize()].is_valid());
                assert!(!binode_map[partner_index.as_usize()].is_valid());
                assert_eq!(
                    node_data,
                    &binode_mapping_function(topology.node_data(partner_index))
                );
                binode_map[node_index.as_usize()] = partner_index.into();
                binode_map[partner_index.as_usize()] = node_index.into();
                data_map.remove(node_data);
            } else {
                let partner_data = binode_mapping_function(node_data);
                assert_ne!(&partner_data, node_data);
                assert_eq!(None, data_map.insert(partner_data, node_index));
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
    <Self as GraphBase>::NodeData: Eq + Hash + Debug,
{
    type Topology = Topology;

    fn new(
        topology: Self::Topology,
        binode_mapping_function: fn(
            &<Self as GraphBase>::NodeData,
        ) -> <Self as GraphBase>::NodeData,
    ) -> Self {
        Self::new_internal(topology, binode_mapping_function, true)
    }

    fn new_unchecked(
        _topology: Self::Topology,
        _binode_mapping_function: fn(
            &<Self as GraphBase>::NodeData,
        ) -> <Self as GraphBase>::NodeData,
    ) -> Self {
        Self::new_internal(_topology, _binode_mapping_function, false)
    }
}

impl<Topology: DynamicGraph> DynamicBigraph for NodeBigraphWrapper<Topology>
where
    <Self as GraphBase>::NodeData: BidirectedData,
    <Self as GraphBase>::EdgeData: Clone,
{
    fn add_partner_nodes(&mut self) {
        for node_id in self.node_indices() {
            if self.partner_node(node_id).is_none() {
                let partner_index = self.add_node(self.node_data(node_id).reverse_complement());
                self.binode_map[node_id.as_usize()] = partner_index.into();
                self.binode_map[partner_index.as_usize()] = node_id.into();
            }
        }
    }

    fn set_partner_nodes(&mut self, a: Self::NodeIndex, b: Self::NodeIndex) {
        assert_ne!(a, b);
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

    fn contains_edge(&self, from: Self::NodeIndex, to: Self::NodeIndex) -> bool {
        self.topology.contains_edge(from, to)
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
}

impl<'a, Topology: NavigableGraph<'a>> NavigableGraph<'a> for NodeBigraphWrapper<Topology> {
    type OutNeighbors = <Topology as NavigableGraph<'a>>::OutNeighbors;
    type InNeighbors = <Topology as NavigableGraph<'a>>::InNeighbors;

    fn out_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::OutNeighbors {
        self.topology.out_neighbors(node_id)
    }

    fn in_neighbors(&'a self, node_id: Self::NodeIndex) -> Self::InNeighbors {
        self.topology.in_neighbors(node_id)
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

#[cfg(test)]
mod tests {
    use super::NodeBigraphWrapper;
    use crate::interface::{
        dynamic_bigraph::DynamicBigraph, static_bigraph::StaticBigraph,
        static_bigraph::StaticBigraphFromDigraph, BidirectedData,
    };
    use crate::traitgraph::implementation::petgraph_impl;
    use crate::traitgraph::interface::{ImmutableGraphContainer, MutableGraphContainer};

    #[test]
    fn test_bigraph_creation() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });

        assert_eq!(Some(n2), bigraph.partner_node(n1));
        assert_eq!(Some(n1), bigraph.partner_node(n2));
        assert_eq!(Some(n4), bigraph.partner_node(n3));
        assert_eq!(Some(n3), bigraph.partner_node(n4));
    }

    #[test]
    #[should_panic]
    fn test_bigraph_creation_unmapped_node() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_node(4);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
    }

    #[test]
    #[should_panic]
    fn test_bigraph_creation_wrongly_mapped_node() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_node(4);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new(graph, |n| {
            if *n == 4 {
                3
            } else if n % 2 == 0 {
                n + 1
            } else {
                n - 1
            }
        });
    }

    #[test]
    #[should_panic]
    fn test_bigraph_creation_self_mapped_node_without_partner() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_node(4);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new(graph, |n| {
            if *n == 4 {
                4
            } else if n % 2 == 0 {
                n + 1
            } else {
                n - 1
            }
        });
    }

    #[test]
    #[should_panic]
    fn test_bigraph_creation_self_mapped_node_with_partner() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_node(4);
        graph.add_node(5);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new(graph, |n| {
            if *n == 4 {
                4
            } else if n % 2 == 0 {
                n + 1
            } else {
                n - 1
            }
        });
    }

    #[test]
    fn test_bigraph_unchecked_creation() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        let bigraph =
            NodeBigraphWrapper::new_unchecked(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });

        assert_eq!(Some(n2), bigraph.partner_node(n1));
        assert_eq!(Some(n1), bigraph.partner_node(n2));
        assert_eq!(Some(n4), bigraph.partner_node(n3));
        assert_eq!(Some(n3), bigraph.partner_node(n4));
    }

    #[test]
    fn test_bigraph_unchecked_creation_unmapped_node() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        let n5 = graph.add_node(4);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        let bigraph =
            NodeBigraphWrapper::new_unchecked(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });

        assert_eq!(Some(n2), bigraph.partner_node(n1));
        assert_eq!(Some(n1), bigraph.partner_node(n2));
        assert_eq!(Some(n4), bigraph.partner_node(n3));
        assert_eq!(Some(n3), bigraph.partner_node(n4));
        assert_eq!(None, bigraph.partner_node(n5));
    }

    #[test]
    #[should_panic]
    fn test_bigraph_unchecked_creation_wrongly_mapped_node_at_end() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_node(4);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new_unchecked(graph, |n| {
            if *n == 4 {
                3
            } else if n % 2 == 0 {
                n + 1
            } else {
                n - 1
            }
        });
    }

    #[test]
    #[should_panic]
    fn test_bigraph_unchecked_creation_wrongly_mapped_node_at_beginning() {
        let mut graph = petgraph_impl::new();
        graph.add_node(4);
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new_unchecked(graph, |n| {
            if *n == 4 {
                3
            } else if n % 2 == 0 {
                n + 1
            } else {
                n - 1
            }
        });
    }

    #[test]
    #[should_panic]
    fn test_bigraph_unchecked_creation_self_mapped_node_without_partner() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_node(4);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new_unchecked(graph, |n| {
            if *n == 4 {
                4
            } else if n % 2 == 0 {
                n + 1
            } else {
                n - 1
            }
        });
    }

    #[test]
    #[should_panic]
    fn test_bigraph_unchecked_creation_self_mapped_node_with_partner() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_node(4);
        graph.add_node(5);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        NodeBigraphWrapper::new_unchecked(graph, |n| {
            if *n == 4 {
                4
            } else if n % 2 == 0 {
                n + 1
            } else {
                n - 1
            }
        });
    }

    #[test]
    fn test_bigraph_verification() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        assert!(bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_verification_self_mapped_node() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        let mut bigraph =
            NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        bigraph.topology.add_node(4);
        bigraph.binode_map.push(4usize.into());
        assert!(!bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_verification_self_unmapped_node() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        let mut bigraph =
            NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        bigraph.topology.add_node(4);
        bigraph.binode_map.push(None.into());
        assert!(!bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_verification_wrongly_mapped_node() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        let mut bigraph =
            NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        bigraph.topology.add_node(4);
        bigraph.binode_map.push(3usize.into());
        assert!(!bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_add_partner_nodes() {
        let mut graph = petgraph_impl::new();
        #[derive(Eq, PartialEq, Debug, Hash, Clone)]
        struct NodeData(u32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
                Self(1000 - self.0)
            }
        }
        let n0 = graph.add_node(NodeData(0));
        let n1 = graph.add_node(NodeData(1));
        graph.add_node(NodeData(2));
        graph.add_node(NodeData(3));
        graph.add_node(NodeData(997));
        graph.add_edge(n0, n1, ());
        let mut graph = NodeBigraphWrapper::new_unchecked(graph, NodeData::reverse_complement);
        assert!(!graph.verify_node_pairing());
        graph.add_partner_nodes();
        assert!(graph.verify_node_pairing());
        assert_eq!(graph.node_count(), 8);
    }
}
