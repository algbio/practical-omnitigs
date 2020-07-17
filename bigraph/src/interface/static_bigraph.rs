use crate::{NodeIndex, StaticGraph};
use num_traits::PrimInt;

/**
 * A node-centric bidirected graph.
 * That is a graph in which each node has a unique partner, and this relation is symmetric.
 */
pub trait StaticBigraph<NodeData, EdgeData, IndexType: PrimInt>:
    StaticGraph<NodeData, EdgeData, IndexType> + Sized
{
    /**
     * Returns the unique partner of the given node id, or `None` if the given node id does not exist.
     */
    fn partner_node(&self, node_id: NodeIndex<IndexType>) -> Option<NodeIndex<IndexType>>;

    /**
     * Returns true if each node has exactly one partner, and this relation is symmetric.
     */
    fn verify_node_pairing(&self) -> bool {
        for node_index in self.node_indices() {
            let partner_index = if let Some(partner_node) = self.partner_node(node_index) {
                partner_node
            } else {
                return false;
            };
            let partner_partner_index =
                if let Some(partner_partner_node) = self.partner_node(partner_index) {
                    partner_partner_node
                } else {
                    return false;
                };

            assert!(
                !node_index.is_invalid()
                    && !partner_index.is_invalid()
                    && !partner_partner_index.is_invalid()
            );
            if node_index != partner_partner_index || node_index == partner_index {
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
    fn verify_mirror_property(&self) -> bool {
        for from_node in self.node_indices() {
            for to_node in self.out_neighbors(from_node).unwrap() {
                let from_node_partner = self.partner_node(from_node).unwrap();
                let to_node_partner = self.partner_node(to_node.node_id).unwrap();
                if !self.contains_edge(to_node_partner, from_node_partner) {
                    return false;
                }
            }
        }

        true
    }

    /**
     * Returns true if all unitigs in the graph have length of at most one node.
     */
    fn verify_unitig_length_is_zero(&self) -> bool {
        for node in self.node_indices() {
            if self.out_neighbors(node).unwrap().into_iter().count() == 1
                && self.in_neighbors(node).unwrap().into_iter().count() == 1
            {
                return false;
            }
        }

        true
    }
}

/**
 * A static bigraph that can be created from a static digraph.
 * Since the graph is static, the resulting topology will be the input topology, only the
 * bigraph node mapping function will be computed on top.
 */
pub trait StaticBigraphFromDigraph<NodeData, EdgeData, IndexType: PrimInt>:
    StaticBigraph<NodeData, EdgeData, IndexType> + Sized
{
    /** The type of directed topology the bigraph is created from. */
    type Topology: StaticGraph<NodeData, EdgeData, IndexType>;

    /**
     * Converts the given topology into a bigraph with the given mapping function.
     * If the resulting graph has wrongly mapped nodes, the method panics.
     */
    fn new(_topology: Self::Topology, _binode_mapping_function: fn(&NodeData) -> NodeData) -> Self;

    /**
     * Converts the given topology into a bigraph with the given mapping function.
     * Wrongly mapped nodes are stored without mapping.
     */
    fn new_unchecked(
        _topology: Self::Topology,
        _binode_mapping_function: fn(&NodeData) -> NodeData,
    ) -> Self;
}

#[cfg(test)]
mod test {
    use crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
    use crate::interface::static_bigraph::StaticBigraph;
    use crate::interface::MutableGraphContainer;
    use crate::{petgraph_impl, StaticBigraphFromDigraph};

    #[test]
    fn test_verify_mirror_property_positive() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        graph.add_edge(n1, n3, 10);
        graph.add_edge(n4, n2, 11);
        graph.add_edge(n3, n1, 12);
        graph.add_edge(n2, n4, 13);
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_mirror_property());
    }

    #[test]
    fn test_verify_mirror_property_unpaired() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        graph.add_edge(n1, n3, 10);
        graph.add_edge(n4, n2, 11);
        graph.add_edge(n3, n1, 12);
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        assert!(bigraph.verify_node_pairing());
        assert!(!bigraph.verify_mirror_property());
    }

    #[test]
    fn test_verify_unitig_length_is_zero_positive() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        let n5 = graph.add_node(4);
        let n6 = graph.add_node(5);
        let n7 = graph.add_node(6);
        let n8 = graph.add_node(7);
        graph.add_edge(n1, n3, 10);
        graph.add_edge(n4, n2, 11);
        graph.add_edge(n3, n5, 12);
        graph.add_edge(n6, n4, 13);
        graph.add_edge(n5, n7, 14);
        graph.add_edge(n8, n6, 15);
        graph.add_edge(n7, n1, 16);
        graph.add_edge(n2, n8, 17);
        graph.add_edge(n1, n5, 18);
        graph.add_edge(n6, n2, 19);
        graph.add_edge(n3, n7, 20);
        graph.add_edge(n8, n4, 21);
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_mirror_property());
        assert!(bigraph.verify_unitig_length_is_zero());
    }

    #[test]
    fn test_verify_unitig_length_is_zero_negative() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        let n5 = graph.add_node(4);
        let n6 = graph.add_node(5);
        let n7 = graph.add_node(6);
        let n8 = graph.add_node(7);
        graph.add_edge(n1, n3, 10);
        graph.add_edge(n4, n2, 11);
        graph.add_edge(n3, n5, 12);
        graph.add_edge(n6, n4, 13);
        graph.add_edge(n5, n7, 14);
        graph.add_edge(n8, n6, 15);
        graph.add_edge(n7, n1, 16);
        graph.add_edge(n2, n8, 17);
        graph.add_edge(n1, n5, 18);
        graph.add_edge(n6, n2, 19);
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_mirror_property());
        assert!(!bigraph.verify_unitig_length_is_zero());
    }
}
