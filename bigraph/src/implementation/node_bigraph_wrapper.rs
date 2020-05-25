use crate::{StaticGraph, EdgeIndex, NodeBigraph, ImmutableGraphContainer, NodeIndices, EdgeIndices};
use crate::NodeIndex;
use std::marker::PhantomData;
use std::hash::Hash;
use std::collections::HashMap;
use num_traits::{PrimInt, NumCast};
use std::fmt::Debug;

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
pub struct NodeBigraphWrapper<NodeData, EdgeData, IndexType: PrimInt, T: StaticGraph<NodeData, EdgeData, IndexType>> {
    pub topology: T,
    binode_map: Vec<NodeIndex<IndexType>>,
    // biedge_map: Vec<EdgeIndex<IndexType>>,
    _p1: PhantomData<NodeData>,
    _p2: PhantomData<EdgeData>,
}

impl<NodeData: Hash + Eq + Debug, EdgeData, IndexType: PrimInt + Debug, T: StaticGraph<NodeData, EdgeData, IndexType>> NodeBigraphWrapper<NodeData, EdgeData, IndexType, T> {
    pub fn new(topology: T, binode_mapping_function: fn(&NodeData) -> NodeData) -> Self {
        let mut data_map: HashMap<NodeData, NodeIndex<IndexType>> = HashMap::new();
        let mut binode_map = vec![NodeIndex::<IndexType>::invalid(); topology.node_count()];
        // let mut biedge_map = vec![NodeIndex::<IndexType>::invalid(); topology.node_count()];

        for node_index in topology.node_indices() {
            let node_data = topology.node_data(node_index).unwrap();

            if let Some(partner_index) = data_map.get(node_data).cloned() {
                assert_eq!(NodeIndex::<IndexType>::invalid(), binode_map[node_index]);
                assert_eq!(NodeIndex::<IndexType>::invalid(), binode_map[partner_index]);
                assert_eq!(node_data, &binode_mapping_function(topology.node_data(partner_index).unwrap()));
                binode_map[node_index] = partner_index;
                binode_map[partner_index] = node_index;
                data_map.remove(node_data);
            } else {
                let partner_data = binode_mapping_function(node_data);
                data_map.insert(partner_data, node_index);
            }
        }

        Self {topology, binode_map, _p1: Default::default(), _p2: Default::default()}

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
}

impl<NodeData, EdgeData, IndexType: PrimInt, T: StaticGraph<NodeData, EdgeData, IndexType>> NodeBigraph<NodeData, EdgeData, IndexType> for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T> {
    fn reverse_complement_node(&self, node_id: NodeIndex<IndexType>) -> Option<NodeIndex<IndexType>> {
        self.binode_map.get(<usize as NumCast>::from(node_id).unwrap()).cloned()
    }
}

impl<NodeData, EdgeData, IndexType: PrimInt, T: StaticGraph<NodeData, EdgeData, IndexType>> ImmutableGraphContainer<NodeData, EdgeData, IndexType> for NodeBigraphWrapper<NodeData, EdgeData, IndexType, T> {
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
}

#[cfg(test)]
mod tests {
    use crate::{petgraph_impl, NodeBigraph, MutableGraphContainer};
    use crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper;

    #[test]
    fn test_bigraph_creation() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        graph.add_edge(n1, n2, "e1"); // Just to fix the EdgeData type parameter
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 {n + 1} else {n - 1});

        assert_eq!(Some(n2), bigraph.reverse_complement_node(n1));
        assert_eq!(Some(n1), bigraph.reverse_complement_node(n2));
        assert_eq!(Some(n4), bigraph.reverse_complement_node(n3));
        assert_eq!(Some(n3), bigraph.reverse_complement_node(n4));
    }
}