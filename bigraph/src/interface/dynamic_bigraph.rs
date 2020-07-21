use crate::interface::static_bigraph::StaticBigraph;
use crate::{DynamicGraph, StaticBigraphFromDigraph};
use num_traits::PrimInt;

pub trait DynamicBigraph<NodeData, EdgeData: Clone, IndexType: PrimInt>:
    DynamicGraph<NodeData, EdgeData, IndexType> + StaticBigraph<NodeData, EdgeData, IndexType>
{
    /**
     * Adds edges such that the graph fulfils the [mirror property].
     * Mirror edges get the cloned `EdgeData` from their existing mirror.
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn add_mirror_edges(&mut self) {
        let mut edges = Vec::new();
        for from_id in self.node_indices() {
            for neighbor in self.out_neighbors(from_id).unwrap() {
                let to_id = neighbor.node_id;
                let mirror_from_id = self.partner_node(to_id).unwrap();
                let mirror_to_id = self.partner_node(from_id).unwrap();
                if !self.contains_edge(mirror_from_id, mirror_to_id) {
                    edges.push((
                        mirror_from_id,
                        self.edge_data(neighbor.edge_id).unwrap().clone(),
                        mirror_to_id,
                    ));
                }
            }
        }

        for edge in edges {
            self.add_edge(edge.0, edge.2, edge.1);
        }
    }
    /**
     * Adds nodes such that the graph becomes a valid bigraph.
     * The indices of existing nodes are not altered.
     */
    fn add_partner_nodes(&mut self);
}

pub trait DynamicBigraphFromDigraph<NodeData, EdgeData: Clone, IndexType: PrimInt>:
    DynamicBigraph<NodeData, EdgeData, IndexType>
    + StaticBigraphFromDigraph<NodeData, EdgeData, IndexType>
where
    Self::Topology: DynamicGraph<NodeData, EdgeData, IndexType>,
{
    /**
     * Converts the given topology into a bigraph.
     * Missing partners are added.
     */
    fn new_with_completed_nodes(
        _topology: Self::Topology,
        _binode_mapping_function: fn(&NodeData) -> NodeData,
    ) -> Self {
        unimplemented!()
    }

    /**
     * Converts the given topology into a bigraph that fulfils the [mirror property].
     * Missing partners and mirror edges are added.
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn new_with_completed_nodes_and_mirrored_edges(
        topology: Self::Topology,
        binode_mapping_function: fn(&NodeData) -> NodeData,
    ) -> Self {
        let mut bigraph = Self::new_with_completed_nodes(topology, binode_mapping_function);
        bigraph.add_mirror_edges();
        bigraph
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        petgraph_impl, BidirectedNodeData, DynamicBigraph, ImmutableGraphContainer,
        MutableGraphContainer, NodeBigraphWrapper, StaticBigraph, StaticBigraphFromDigraph,
    };

    #[test]
    fn test_bigraph_add_mirror_edges() {
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
        let n2 = graph.add_node(NodeData(2));
        let n3 = graph.add_node(NodeData(3));
        let n4 = graph.add_node(NodeData(997));
        graph.add_edge(n3, n4, ()); // This edge is a self-mirror
        graph.add_edge(n1, n2, ()); // This edge is not a self-mirror
        graph.add_edge(n0, n3, ()); // This edge is not a self-mirror

        let mut graph = NodeBigraphWrapper::new_unchecked(graph, NodeData::reverse_complement);
        graph.add_partner_nodes();
        assert!(graph.verify_node_pairing());
        assert_eq!(graph.node_count(), 8);

        assert!(!graph.verify_mirror_property());
        graph.add_mirror_edges();
        assert!(graph.verify_mirror_property());
        assert_eq!(graph.edge_count(), 5);
    }
}
