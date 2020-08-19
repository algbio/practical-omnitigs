use crate::interface::BidirectedData;
use traitgraph::interface::{GraphBase, StaticGraph};

/**
 * A node-centric bidirected graph.
 * That is a graph in which each node has a unique partner, and this relation is symmetric.
 */
pub trait StaticBigraph: StaticGraph {
    /**
     * Returns the unique partner of the given node id, or `None` if the given node id has no partner node.
     */
    fn partner_node(&self, node_id: Self::NodeIndex) -> Option<Self::NodeIndex>;

    /**
     * Returns true if each node has exactly one partner, and this relation is symmetric.
     * This check allows nodes that are their own partner.
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
            if node_index != partner_partner_index {
                return false;
            }
        }

        true
    }

    /**
     * Returns true if each node has exactly one partner, and this relation is symmetric and irreflexive (no node is its own partner).
     */
    fn verify_node_pairing_without_self_partners(&self) -> bool {
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
            if node_index != partner_partner_index || node_index == partner_index {
                return false;
            }
        }

        true
    }
}

/**
 * A edge-centric bidirected graph.
 * That is a graph in which each node has a unique partner, and this relation is symmetric.
 */
pub trait StaticNodeCentricBigraph: StaticBigraph {
    /**
     * Returns the unique partner of the given edge id, or `None` if the given edge id has no partner edge.
     * If the edge is its own reverse complement, and an partner edge with a different id exists, then the different id is returned.
     * Otherwise, for an edge that is its own reverse complement, the given id is returned.
     */
    fn partner_edge_node_centric(&self, edge_id: Self::EdgeIndex) -> Option<Self::EdgeIndex> {
        let endpoints = self.edge_endpoints(edge_id);
        let reverse_from = self.partner_node(endpoints.to_node)?;
        let reverse_to = self.partner_node(endpoints.from_node)?;
        let mut result = None;

        for reverse_edge in self.out_neighbors_to(reverse_from, reverse_to) {
            debug_assert_eq!(reverse_edge.node_id, reverse_to);
            let reverse_edge_id = reverse_edge.edge_id;
            if let Some(node) = result {
                if node == edge_id {
                    return Some(reverse_edge_id);
                }
            } else if reverse_edge_id == edge_id {
                result = Some(reverse_edge_id);
            } else {
                return Some(reverse_edge_id);
            }
        }

        None
    }

    /**
     * Returns true if the node-centric [mirror property] of edges is fulfilled.
     * Assumes that the node pairing is correct (See [verify_node_pairing()](NodeBigraphWrapper::verify_node_pairing))
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn verify_node_mirror_property(&self) -> bool {
        for from_node in self.node_indices() {
            for to_node in self.out_neighbors(from_node) {
                let from_node_partner = self.partner_node(from_node).unwrap();
                let to_node_partner = self.partner_node(to_node.node_id).unwrap();
                if self.edge_count_between(to_node_partner, from_node_partner)
                    != self.edge_count_between(from_node, to_node.node_id)
                {
                    return false;
                }
            }
        }

        true
    }
}

/**
 * A edge-centric bidirected graph.
 * That is a graph in which each node and each edge has a unique partner, and this relation is symmetric.
 */
pub trait StaticEdgeCentricBigraph: StaticBigraph
where
    <Self as GraphBase>::EdgeData: BidirectedData + Eq,
{
    /**
     * Returns the unique partner of the given edge id, or `None` if the given edge id has no partner edge.
     * If the edge is its own reverse complement, and an partner edge with a different id exists, then the different id is returned.
     * Otherwise, for an edge that is its own reverse complement, the given id is returned.
     */
    fn partner_edge_edge_centric(&self, edge_id: Self::EdgeIndex) -> Option<Self::EdgeIndex> {
        let endpoints = self.edge_endpoints(edge_id);
        let reverse_from = self.partner_node(endpoints.to_node)?;
        let reverse_to = self.partner_node(endpoints.from_node)?;
        let edge_data = self.edge_data(edge_id);
        let mut result = None;

        for reverse_edge in self.out_neighbors_to(reverse_from, reverse_to) {
            debug_assert_eq!(reverse_edge.node_id, reverse_to);
            let reverse_edge_id = reverse_edge.edge_id;
            if &edge_data.reverse_complement() == self.edge_data(reverse_edge_id) {
                if let Some(node) = result {
                    if node == edge_id {
                        return Some(reverse_edge_id);
                    }
                } else if reverse_edge_id == edge_id {
                    result = Some(reverse_edge_id);
                } else {
                    return Some(reverse_edge_id);
                }
            }
        }

        None
    }

    /**
     * Returns true if the edge-centric [mirror property] of edges is fulfilled.
     * Assumes that the node pairing is correct (See [verify_node_pairing()](NodeBigraphWrapper::verify_node_pairing))
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn verify_edge_mirror_property(&self) -> bool {
        for from_node in self.node_indices() {
            for neighbor in self.out_neighbors(from_node) {
                let edge = neighbor.edge_id;
                if self.partner_edge_edge_centric(edge).is_none() {
                    return false;
                }
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
pub trait StaticBigraphFromDigraph: StaticBigraph {
    /** The type of directed topology the bigraph is created from. */
    type Topology: StaticGraph<NodeData = Self::NodeData, EdgeData = Self::EdgeData>;

    /**
     * Converts the given topology into a bigraph with the given mapping function.
     * If the resulting graph has wrongly mapped nodes, the method panics.
     */
    fn new(
        _topology: Self::Topology,
        _binode_mapping_function: fn(&Self::NodeData) -> Self::NodeData,
    ) -> Self;

    /**
     * Converts the given topology into a bigraph with the given mapping function.
     * Unmapped nodes are stored without mapping.
     * If a node maps to another node, but the other node does not map back, then this method panics.
     */
    fn new_unchecked(
        _topology: Self::Topology,
        _binode_mapping_function: fn(&Self::NodeData) -> Self::NodeData,
    ) -> Self;
}

#[cfg(test)]
mod test {
    use crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
    use crate::interface::static_bigraph::StaticBigraph;
    use crate::interface::static_bigraph::StaticBigraphFromDigraph;
    use crate::interface::static_bigraph::StaticNodeCentricBigraph;
    use crate::interface::BidirectedData;
    use crate::traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::MutableGraphContainer;

    #[derive(Debug, Clone, Copy, Eq, PartialEq)]
    struct EdgeData(usize);

    impl BidirectedData for EdgeData {
        fn reverse_complement(&self) -> Self {
            EdgeData(1000 - self.0)
        }
    }

    #[test]
    fn test_verify_node_mirror_property_positive() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(11));
        graph.add_edge(n3, n1, EdgeData(12));
        graph.add_edge(n2, n4, EdgeData(13));
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_node_mirror_property());
    }

    #[test]
    fn test_verify_node_mirror_property_unpaired() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(11));
        graph.add_edge(n3, n1, EdgeData(12));
        let bigraph = NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        assert!(bigraph.verify_node_pairing());
        assert!(!bigraph.verify_node_mirror_property());
    }

    #[test]
    fn test_verify_node_mirror_property_duplicate_edges() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(0);
        let n2 = graph.add_node(1);
        let n3 = graph.add_node(2);
        let n4 = graph.add_node(3);
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(11));
        graph.add_edge(n3, n1, EdgeData(12));
        graph.add_edge(n2, n4, EdgeData(13));

        let mut bigraph =
            NodeBigraphWrapper::new(graph, |n| if n % 2 == 0 { n + 1 } else { n - 1 });
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_node_mirror_property());

        bigraph.add_edge(n1, n3, EdgeData(14));
        assert!(bigraph.verify_node_pairing());
        assert!(!bigraph.verify_node_mirror_property());

        bigraph.add_edge(n4, n2, EdgeData(15));
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_node_mirror_property());
    }
}
