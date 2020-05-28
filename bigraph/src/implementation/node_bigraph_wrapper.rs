use crate::{NodeIndex, NavigableGraph};
use crate::{
    EdgeIndex, EdgeIndices, ImmutableGraphContainer, NodeBigraph, NodeIndices, StaticGraph,
};
use num_traits::{NumCast, PrimInt, ToPrimitive};
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::marker::PhantomData;

/**
*   Wrapper for a static graph that adds a binode mapping function.
*
*   Bigraphs can be represented with this struct by creating their topology as normal directed graph where each binode is split into its two parts.
*   The binode mapping function then associates the parts with each other.
*
*   ```rust
*   use bigraph::node_bigraph_wrapper::NodeBigraphWrapper;
*   use bigraph::{NodeBigraph, MutableGraphContainer};
*   use bigraph::petgraph_impl;
*
*   let mut graph = petgraph_impl::new();
*   let n1 = graph.add_node(0);
*   let n2 = graph.add_node(1);
*   graph.add_edge(n1.clone(), n2.clone(), ());
*   graph.add_edge(n2.clone(), n1.clone(), ());
*   let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 {n + 1} else {n - 1});
*   assert_eq!(Some(n2.clone()), bigraph.reverse_complement_node(n1.clone()));
*   assert_eq!(Some(n1.clone()), bigraph.reverse_complement_node(n2.clone()));
*   ```
*/
pub struct NodeBigraphWrapper<
    NodeData,
    EdgeData,
    IndexType: PrimInt,
    T,
> {
    pub topology: T,
    binode_map: Vec<NodeIndex<IndexType>>,
    // biedge_map: Vec<EdgeIndex<IndexType>>,
    _p1: PhantomData<NodeData>,
    _p2: PhantomData<EdgeData>,
}

impl<
        NodeData: Hash + Eq + Debug,
        EdgeData,
        IndexType: PrimInt + Debug,
        T: StaticGraph<NodeData, EdgeData, IndexType>,
    > NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
    pub fn new(topology: T, binode_mapping_function: fn(&NodeData) -> NodeData) -> Self {
        let mut data_map: HashMap<NodeData, NodeIndex<IndexType>> = HashMap::new();
        let mut binode_map = vec![NodeIndex::<IndexType>::invalid(); topology.node_count()];
        // let mut biedge_map = vec![NodeIndex::<IndexType>::invalid(); topology.node_count()];

        for node_index in topology.node_indices() {
            let node_data = topology.node_data(node_index).unwrap();

            if let Some(partner_index) = data_map.get(node_data).cloned() {
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
                data_map.insert(partner_data, node_index);
            }
        }

        assert!(data_map.is_empty());
        Self {
            topology,
            binode_map,
            _p1: Default::default(),
            _p2: Default::default(),
        }

        /*for edge_index in topology.edge_indices() {
            let edge = topology.edge(edge_index);

            // Search reverse edge
            let partner_from = binode_mapping_function(topology.node_data(edge.to).unwrap());
            let partner_to = binode_mapping_function(topology.node_data(edge.from).unwrap());

            for neighbor_node in topology.out_neighbors(partner_from) {

            }
        }
        unimplemented!()
        */
    }

    /**
     * Returns true if each node has exactly one partner, and this relation is symmetric.
     */
    pub fn verify_node_pairing(&self) -> bool {
        if self.binode_map.len() != self.topology.node_count() {
            return false;
        }

        for (node_index, partner_index) in self.binode_map.iter().enumerate() {
            let node_index_inner = if let Some(node_index_inner) = IndexType::from(node_index) {
                node_index_inner
            } else {
                return false;
            };
            let node_index = NodeIndex::from(node_index_inner);
            let partner_index_inner = if let Some(parner_index_inner) = partner_index.to_usize() {
                parner_index_inner
            } else {
                return false;
            };
            let reverse_mapping =
                if let Some(reverse_mapping) = self.binode_map.get(partner_index_inner).cloned() {
                    reverse_mapping
                } else {
                    return false;
                };
            if reverse_mapping != node_index
                || node_index == *partner_index
                || reverse_mapping.is_invalid()
                || partner_index.is_invalid()
                || node_index.is_invalid()
            {
                return false;
            }
        }

        true
    }

    /**
     * Returns true if the [mirror property] of edges is fulfilled.
     * Assumes that the node pairing is correct (See [verify_node_pairing()](NodeBigraphWrapper::verify_node_pairing))
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    pub fn verify_mirror_property(&self) -> bool {
        for from_node in self.node_indices() {
            for to_node in self.out_neighbors(from_node) {
                let from_node_partner = self.partner_node(from_node).unwrap();
                let to_node_partner = self.partner_node(to_node).unwrap();
                if !self.contains_edge(to_node_partner, from_node_partner) {
                    return false
                }
            }
        }

        true
    }

    /**
     * Returns true if all unitigs in the graph have length of at most one node.
     */
    pub fn verify_unitig_length_is_zero(&self) -> bool {
        unimplemented!()
    }
}

impl<NodeData, EdgeData, IndexType: PrimInt, T>
    NodeBigraph<NodeData, EdgeData, IndexType>
    for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
    fn partner_node(
        &self,
        node_id: NodeIndex<IndexType>,
    ) -> Option<NodeIndex<IndexType>> {
        self.binode_map
            .get(<usize as NumCast>::from(node_id).unwrap())
            .cloned()
    }
}

impl<NodeData, EdgeData, IndexType: PrimInt, T: ImmutableGraphContainer<NodeData, EdgeData, IndexType>>
    ImmutableGraphContainer<NodeData, EdgeData, IndexType>
    for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T>
{
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

impl<'a, NodeData, EdgeData, IndexType: PrimInt, T: NavigableGraph<'a, NodeData, EdgeData, IndexType>>
NavigableGraph<'a, NodeData, EdgeData, IndexType>
for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T> {
    type OutNeighbors = <T as NavigableGraph<'a, NodeData, EdgeData, IndexType>>::OutNeighbors;
    type InNeighbors = <T as NavigableGraph<'a, NodeData, EdgeData, IndexType>>::InNeighbors;

    fn out_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::OutNeighbors> {
        self.topology.out_neighbors(node_id)
    }

    fn in_neighbors(&'a self, node_id: NodeIndex<IndexType>) -> Option<Self::InNeighbors> {
        self.topology.in_neighbors(node_id)
    }
}

#[cfg(test)]
mod tests {
    use crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
    use crate::{petgraph_impl, MutableGraphContainer, NodeBigraph, NodeIndex};

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
}
