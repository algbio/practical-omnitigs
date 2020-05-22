use crate::StaticGraph;
use crate::NodeIndex;
use std::marker::PhantomData;
use std::hash::Hash;
use std::collections::HashMap;
use num_traits::PrimInt;
use petgraph::dot::Config::NodeIndexLabel;
use std::fmt::Debug;

pub struct BigraphWrapper<NodeData, EdgeData, IndexType: PrimInt, T: StaticGraph<NodeData, EdgeData, IndexType>> {
    topology: T,
    binode_map: Vec<NodeIndex<IndexType>>,
    _p1: PhantomData<NodeData>,
    _p2: PhantomData<EdgeData>,
}

impl<NodeData: Hash + Eq, EdgeData, IndexType: PrimInt + Debug, T: StaticGraph<NodeData, EdgeData, IndexType>> BigraphWrapper<NodeData, EdgeData, IndexType, T> {
    fn new(topology: T, binode_mapping_function: fn(&NodeData) -> NodeData) -> Self {
        let mut data_map: HashMap<NodeData, NodeIndex<IndexType>> = HashMap::new();
        let mut binode_map = vec![NodeIndex::<IndexType>::invalid(); topology.node_count().to_usize().unwrap()];

        for node_index in topology.node_indices() {
            let node_data = topology.node_data(node_index).unwrap();

            if let Some(partner_index) = data_map.get(node_data).cloned() {
                assert_eq!(NodeIndex::<IndexType>::invalid(), binode_map[node_index]);
                assert_eq!(NodeIndex::<IndexType>::invalid(), binode_map[partner_index]);
                binode_map[node_index] = partner_index;
                binode_map[partner_index] = node_index;
                data_map.remove(node_data);
            } else {
                data_map.insert(binode_mapping_function(node_data), node_index);
            }
        }

        unimplemented!()
    }
}