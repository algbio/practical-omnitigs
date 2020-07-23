use crate::{BidirectedNodeData, DynamicBigraph, StaticBigraph, StaticBigraphFromDigraph};
use num_traits::{NumCast, PrimInt};
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::marker::PhantomData;
use traitgraph::{
    DynamicGraph, EdgeIndex, EdgeIndices, ImmutableGraphContainer, MutableGraphContainer,
    NavigableGraph, NodeIndex, NodeIndices, StaticGraph,
};

/**
*   Wrapper for a static graph that adds a binode mapping function.
*
*   Bigraphs can be represented with this struct by creating their topology as normal directed graph where each binode is split into its two parts.
*   The binode mapping function then associates the parts with each other.
*
*   ```rust
*   use bigraph::node_bigraph_wrapper::NodeBigraphWrapper;
*   use bigraph::{StaticBigraph, MutableGraphContainer, StaticBigraphFromDigraph};
*   use bigraph::petgraph_impl;
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
pub struct NodeBigraphWrapper<NodeData, EdgeData, IndexType: PrimInt, Topology> {
    pub topology: Topology,
    binode_map: Vec<NodeIndex<IndexType>>,
    // biedge_map: Vec<EdgeIndex<IndexType>>,
    _p1: PhantomData<NodeData>,
    _p2: PhantomData<EdgeData>,
}

impl<NodeData, EdgeData, IndexType: PrimInt, T: StaticGraph<NodeData, EdgeData, IndexType>>
    StaticBigraph<NodeData, EdgeData, IndexType>
    for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
    fn partner_node(&self, node_id: NodeIndex<IndexType>) -> Option<NodeIndex<IndexType>> {
        if node_id.is_invalid() {
            return None;
        }

        if let Some(partner_node) = self
            .binode_map
            .get(<usize as NumCast>::from(node_id).unwrap())
        {
            if partner_node.is_invalid() {
                None
            } else {
                Some(*partner_node)
            }
        } else {
            None
        }
    }
}

impl<
        NodeData: Eq + Hash + Debug,
        EdgeData,
        IndexType: PrimInt + Debug,
        Topology: StaticGraph<NodeData, EdgeData, IndexType>,
    > NodeBigraphWrapper<NodeData, EdgeData, IndexType, Topology>
{
    fn new_internal(
        topology: Topology,
        binode_mapping_function: fn(&NodeData) -> NodeData,
        checked: bool,
    ) -> Self {
        let mut data_map: HashMap<NodeData, NodeIndex<IndexType>> = HashMap::new();
        let mut binode_map = vec![NodeIndex::<IndexType>::invalid(); topology.node_count()];
        // let mut biedge_map = vec![NodeIndex::<IndexType>::invalid(); topology.node_count()];

        for node_index in topology.node_indices() {
            let node_data = topology.node_data(node_index).unwrap();

            if let Some(partner_index) = data_map.get(node_data).cloned() {
                assert_ne!(node_index, partner_index);
                assert_eq!(NodeIndex::<IndexType>::invalid(), binode_map[node_index]);
                assert_eq!(NodeIndex::<IndexType>::invalid(), binode_map[partner_index]);
                assert_eq!(
                    node_data,
                    &binode_mapping_function(topology.node_data(partner_index).unwrap())
                );
                binode_map[node_index] = partner_index;
                binode_map[partner_index] = node_index;
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
                let node_data = topology.node_data(node_index).unwrap();
                assert!(!data_map.contains_key(node_data));
            }
        }
        Self {
            topology,
            binode_map,
            _p1: Default::default(),
            _p2: Default::default(),
        }
    }
}

impl<
        NodeData: Eq + Hash + Debug,
        EdgeData,
        IndexType: PrimInt + Debug,
        Topology: StaticGraph<NodeData, EdgeData, IndexType>,
    > StaticBigraphFromDigraph<NodeData, EdgeData, IndexType>
    for NodeBigraphWrapper<NodeData, EdgeData, IndexType, Topology>
{
    type Topology = Topology;

    fn new(topology: Self::Topology, binode_mapping_function: fn(&NodeData) -> NodeData) -> Self {
        Self::new_internal(topology, binode_mapping_function, true)
    }

    fn new_unchecked(
        _topology: Self::Topology,
        _binode_mapping_function: fn(&NodeData) -> NodeData,
    ) -> Self {
        Self::new_internal(_topology, _binode_mapping_function, false)
    }
}

impl<
        NodeData: BidirectedNodeData,
        EdgeData: Clone,
        IndexType: PrimInt,
        T: DynamicGraph<NodeData, EdgeData, IndexType>,
    > DynamicBigraph<NodeData, EdgeData, IndexType>
    for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
    fn add_partner_nodes(&mut self) {
        for node_id in self.node_indices() {
            if self.partner_node(node_id).is_none() {
                let partner_index =
                    self.add_node(self.node_data(node_id).unwrap().reverse_complement());
                self.binode_map[node_id] = partner_index;
                self.binode_map[<usize as NumCast>::from(partner_index).unwrap()] = node_id;
            }
        }
    }
}

impl<
        NodeData,
        EdgeData,
        IndexType: PrimInt,
        T: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>,
    > ImmutableGraphContainer for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
    type NodeData = NodeData;
    type EdgeData = EdgeData;
    type IndexType = IndexType;

    fn node_indices(&self) -> NodeIndices<IndexType> {
        self.topology.node_indices()
    }

    fn edge_indices(&self) -> EdgeIndices<IndexType> {
        self.topology.edge_indices()
    }

    fn node_count(&self) -> usize {
        self.topology.node_count()
    }

    fn edge_count(&self) -> usize {
        self.topology.edge_count()
    }

    fn node_data(&self, node_id: NodeIndex<IndexType>) -> Option<&NodeData> {
        self.topology.node_data(node_id)
    }

    fn edge_data(&self, edge_id: EdgeIndex<IndexType>) -> Option<&EdgeData> {
        self.topology.edge_data(edge_id)
    }

    fn node_data_mut(&mut self, node_id: NodeIndex<IndexType>) -> Option<&mut NodeData> {
        self.topology.node_data_mut(node_id)
    }

    fn edge_data_mut(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<&mut EdgeData> {
        self.topology.edge_data_mut(edge_id)
    }

    fn contains_edge(&self, from: NodeIndex<IndexType>, to: NodeIndex<IndexType>) -> bool {
        self.topology.contains_edge(from, to)
    }
}

impl<
        NodeData,
        EdgeData,
        IndexType: PrimInt,
        T: MutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData, IndexType = IndexType>,
    > MutableGraphContainer for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
    type NodeData = NodeData;
    type EdgeData = EdgeData;
    type IndexType = IndexType;

    fn add_node(&mut self, node_data: NodeData) -> NodeIndex<IndexType> {
        self.binode_map.push(NodeIndex::invalid());
        self.topology.add_node(node_data)
    }

    fn add_edge(
        &mut self,
        from: NodeIndex<IndexType>,
        to: NodeIndex<IndexType>,
        edge_data: EdgeData,
    ) -> EdgeIndex<IndexType> {
        self.topology.add_edge(from, to, edge_data)
    }

    fn remove_node(&mut self, node_id: NodeIndex<IndexType>) -> Option<NodeData> {
        self.topology.remove_node(node_id)
    }

    fn remove_edge(&mut self, edge_id: EdgeIndex<IndexType>) -> Option<EdgeData> {
        self.topology.remove_edge(edge_id)
    }
}

impl<
        'a,
        NodeData,
        EdgeData,
        IndexType: PrimInt,
        T: NavigableGraph<'a, NodeData, EdgeData, IndexType>,
    > NavigableGraph<'a, NodeData, EdgeData, IndexType>
    for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
    type OutNeighbors = <T as NavigableGraph<'a, NodeData, EdgeData, IndexType>>::OutNeighbors;
    type InNeighbors = <T as NavigableGraph<'a, NodeData, EdgeData, IndexType>>::InNeighbors;

    fn out_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::OutNeighbors> {
        self.topology.out_neighbors(node_id)
    }

    fn in_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::InNeighbors> {
        self.topology.in_neighbors(node_id)
    }
}

impl<NodeData, EdgeData, IndexType: PrimInt, T: Default> Default
    for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
    fn default() -> Self {
        Self {
            topology: T::default(),
            binode_map: vec![],
            _p1: Default::default(),
            _p2: Default::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
    use crate::{
        petgraph_impl, BidirectedNodeData, DynamicBigraph, ImmutableGraphContainer,
        MutableGraphContainer, NodeIndex, StaticBigraph, StaticBigraphFromDigraph,
    };

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
        bigraph.binode_map.push(NodeIndex::from(4usize));
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
        bigraph.binode_map.push(NodeIndex::invalid());
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
        bigraph.binode_map.push(NodeIndex::from(3usize));
        assert!(!bigraph.verify_node_pairing());
    }

    #[test]
    fn test_bigraph_add_partner_nodes() {
        let mut graph = petgraph_impl::new();
        #[derive(Eq, PartialEq, Debug, Hash, Clone)]
        struct NodeData(u32);
        impl BidirectedNodeData for NodeData {
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
