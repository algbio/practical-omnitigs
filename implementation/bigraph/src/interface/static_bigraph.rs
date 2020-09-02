use crate::interface::BidirectedData;
use std::collections::HashSet;
use traitgraph::interface::{GraphBase, StaticGraph};

/**
 * A node-centric bidirected graph.
 * That is a graph in which each node has a unique mirror, and this relation is symmetric.
 */
pub trait StaticBigraph: StaticGraph {
    /**
     * Returns the unique mirror of the given node id, or `None` if the given node id has no mirror node.
     */
    fn mirror_node(&self, node_id: Self::NodeIndex) -> Option<Self::NodeIndex>;

    /**
     * Returns true if each node has exactly one mirror, and this relation is symmetric.
     * This check allows nodes that are their own mirror.
     */
    fn verify_node_pairing(&self) -> bool {
        for node_index in self.node_indices() {
            let mirror_index = if let Some(mirror_node) = self.mirror_node(node_index) {
                mirror_node
            } else {
                return false;
            };
            let mirror_mirror_index =
                if let Some(mirror_mirror_node) = self.mirror_node(mirror_index) {
                    mirror_mirror_node
                } else {
                    return false;
                };
            if node_index != mirror_mirror_index {
                return false;
            }
        }

        true
    }

    /**
     * Returns true if each node has exactly one mirror, and this relation is symmetric and irreflexive (no node is its own mirror).
     */
    fn verify_node_pairing_without_self_mirrors(&self) -> bool {
        for node_index in self.node_indices() {
            let mirror_index = if let Some(mirror_node) = self.mirror_node(node_index) {
                mirror_node
            } else {
                return false;
            };
            let mirror_mirror_index =
                if let Some(mirror_mirror_node) = self.mirror_node(mirror_index) {
                    mirror_mirror_node
                } else {
                    return false;
                };
            if node_index != mirror_mirror_index || node_index == mirror_index {
                return false;
            }
        }

        true
    }
}

/**
 * A edge-centric bidirected graph.
 * That is a graph in which each node has a unique mirror, and this relation is symmetric.
 */
pub trait StaticNodeCentricBigraph: StaticBigraph {
    /**
     * Returns the unique mirror of the given edge id, or `None` if the given edge id has no mirror edge.
     * If the edge is its own reverse complement, and an mirror edge with a different id exists, then the different id is returned.
     * Otherwise, for an edge that is its own reverse complement, the given id is returned.
     */
    fn mirror_edge_node_centric(&self, edge_id: Self::EdgeIndex) -> Option<Self::EdgeIndex> {
        let endpoints = self.edge_endpoints(edge_id);
        let reverse_from = self.mirror_node(endpoints.to_node)?;
        let reverse_to = self.mirror_node(endpoints.from_node)?;
        let mut result = None;

        for reverse_edge_id in self.out_neighbors_to(reverse_from, reverse_to) {
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
                let from_node_mirror = self.mirror_node(from_node).unwrap();
                let to_node_mirror = self.mirror_node(to_node.node_id).unwrap();
                if self.edge_count_between(to_node_mirror, from_node_mirror)
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
 * That is a graph in which each node and each edge has a unique mirror, and this relation is symmetric.
 */
pub trait StaticEdgeCentricBigraph: StaticBigraph
where
    <Self as GraphBase>::EdgeData: BidirectedData + Eq,
{
    /**
     * Returns the unique mirror of the given edge id, or `None` if the given edge id has no mirror edge.
     * If the edge is its own reverse complement, and an mirror edge with a different id exists, then the different id is returned.
     * Otherwise, for an edge that is its own reverse complement, the given id is returned.
     */
    fn mirror_edge_edge_centric(&self, edge_id: Self::EdgeIndex) -> Option<Self::EdgeIndex> {
        let endpoints = self.edge_endpoints(edge_id);
        let reverse_from = self.mirror_node(endpoints.to_node)?;
        let reverse_to = self.mirror_node(endpoints.from_node)?;
        let edge_data = self.edge_data(edge_id);
        let mut result = None;

        for reverse_edge_id in self.out_neighbors_to(reverse_from, reverse_to) {
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
     * Assumes that the node pairing is correct (See [verify_node_pairing()](NodeBigraphWrapper::verify_node_pairing)) and that no two edges are the same, except for self mirrors.
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn verify_edge_mirror_property(&self) -> bool {
        let mut edge_set = HashSet::new();

        for from_node in self.node_indices() {
            for neighbor in self.out_neighbors(from_node) {
                let edge = neighbor.edge_id;
                if let Some(mirror_edge) = self.mirror_edge_edge_centric(edge) {
                    let to_node = neighbor.node_id;
                    let complete_edge = (from_node, to_node, edge);
                    let mirror_complete_edge = (
                        self.mirror_node(to_node).unwrap(),
                        self.mirror_node(from_node).unwrap(),
                        mirror_edge,
                    );

                    if edge_set.contains(&mirror_complete_edge) {
                        edge_set.remove(&mirror_complete_edge);
                        if &self.edge_data(edge).reverse_complement() != self.edge_data(mirror_edge)
                        {
                            return false;
                        }
                    } else {
                        edge_set.insert(complete_edge);
                    }
                } else {
                    return false;
                }
            }
        }

        edge_set.is_empty()
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
    fn new(topology: Self::Topology) -> Self;

    /**
     * Converts the given topology into a bigraph with the given mapping function.
     * Unmapped nodes are stored without mapping.
     * If a node maps to another node, but the other node does not map back, then this method panics.
     */
    fn new_unchecked(topology: Self::Topology) -> Self;
}

#[cfg(test)]
mod test {
    use crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
    use crate::interface::dynamic_bigraph::DynamicBigraph;
    use crate::interface::static_bigraph::StaticBigraph;
    use crate::interface::static_bigraph::StaticBigraphFromDigraph;
    use crate::interface::static_bigraph::StaticEdgeCentricBigraph;
    use crate::interface::static_bigraph::StaticNodeCentricBigraph;
    use crate::interface::BidirectedData;
    use crate::traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::ImmutableGraphContainer;
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
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
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
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(11));
        graph.add_edge(n3, n1, EdgeData(12));
        graph.add_edge(n2, n4, EdgeData(13));
        let bigraph = NodeBigraphWrapper::new(graph);
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_node_mirror_property());
    }

    #[test]
    fn test_verify_node_mirror_property_unpaired() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
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
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(11));
        graph.add_edge(n3, n1, EdgeData(12));
        let bigraph = NodeBigraphWrapper::new(graph);
        assert!(bigraph.verify_node_pairing());
        assert!(!bigraph.verify_node_mirror_property());
    }

    #[test]
    fn test_verify_node_mirror_property_duplicate_edges() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
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
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(11));
        graph.add_edge(n3, n1, EdgeData(12));
        graph.add_edge(n2, n4, EdgeData(13));

        let mut bigraph = NodeBigraphWrapper::new(graph);
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_node_mirror_property());

        bigraph.add_edge(n1, n3, EdgeData(14));
        assert!(bigraph.verify_node_pairing());
        assert!(!bigraph.verify_node_mirror_property());

        bigraph.add_edge(n4, n2, EdgeData(15));
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_node_mirror_property());
    }

    #[test]
    fn test_verify_edge_mirror_property_positive() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
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
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(990));
        graph.add_edge(n3, n1, EdgeData(12));
        graph.add_edge(n2, n4, EdgeData(988));
        let bigraph = NodeBigraphWrapper::new(graph);
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_edge_mirror_property());
    }

    #[test]
    fn test_verify_edge_mirror_property_unpaired() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
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
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(990));
        graph.add_edge(n3, n1, EdgeData(12));
        graph.add_edge(n4, n2, EdgeData(990));
        let bigraph = NodeBigraphWrapper::new(graph);
        assert!(bigraph.verify_node_pairing());
        assert!(!bigraph.verify_edge_mirror_property());
    }

    #[test]
    fn test_verify_edge_mirror_property_duplicate_edges_with_differing_data() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
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
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n4, n2, EdgeData(990));
        graph.add_edge(n3, n1, EdgeData(12));
        graph.add_edge(n2, n4, EdgeData(988));

        let mut bigraph = NodeBigraphWrapper::new(graph);
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_edge_mirror_property());

        bigraph.add_edge(n1, n3, EdgeData(14));
        assert!(bigraph.verify_node_pairing());
        assert!(!bigraph.verify_edge_mirror_property());

        bigraph.add_edge(n4, n2, EdgeData(986));
        assert!(bigraph.verify_node_pairing());
        assert!(bigraph.verify_edge_mirror_property());
    }

    #[test]
    fn test_verify_edge_mirror_property_duplicate_edges_with_plus_minus_loop() {
        #[derive(Clone, Eq, PartialEq, Hash, Debug)]
        struct NodeData(i32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
                Self(1000 - self.0)
            }
        }

        let mut graph = NodeBigraphWrapper::new(petgraph_impl::new());
        let n1 = graph.add_node(NodeData(0));
        let n2 = graph.add_node(NodeData(1000));
        let n3 = graph.add_node(NodeData(500));
        graph.set_mirror_nodes(n1, n2);
        graph.set_mirror_nodes(n3, n3);
        graph.add_edge(n1, n3, EdgeData(10));
        graph.add_edge(n3, n2, EdgeData(990));
        graph.add_edge(n3, n1, EdgeData(12));
        graph.add_edge(n2, n3, EdgeData(988));

        assert!(graph.verify_node_pairing());
        assert!(graph.verify_edge_mirror_property());

        graph.add_edge(n1, n3, EdgeData(14));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n3, n2, EdgeData(986));
        assert!(graph.verify_node_pairing());
        assert!(graph.verify_edge_mirror_property());

        assert_eq!(graph.edge_count(), 6);
        let mirror_copy = graph.clone();

        graph.add_edge(n1, n3, EdgeData(14));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n3, n2, EdgeData(986));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());
        assert_eq!(graph.edge_count(), 8);

        let mut graph = mirror_copy.clone();

        graph.add_edge(n1, n3, EdgeData(100));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n1, n3, EdgeData(100));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n3, n2, EdgeData(900));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n3, n2, EdgeData(900));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());
        assert_eq!(graph.edge_count(), 10);

        let mut graph = mirror_copy.clone();

        graph.add_edge(n3, n2, EdgeData(900));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n1, n3, EdgeData(100));
        assert!(graph.verify_node_pairing());
        assert!(graph.verify_edge_mirror_property());

        graph.add_edge(n1, n3, EdgeData(100));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n3, n2, EdgeData(900));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());
        assert_eq!(graph.edge_count(), 10);

        let mut graph = mirror_copy.clone();

        graph.add_edge(n3, n2, EdgeData(986));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n1, n3, EdgeData(14));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());
        assert_eq!(graph.edge_count(), 8);

        let mut graph = mirror_copy;

        graph.add_edge(n3, n3, EdgeData(500));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n3, n3, EdgeData(500));
        assert!(graph.verify_node_pairing());
        assert!(graph.verify_edge_mirror_property());
        assert_eq!(graph.edge_count(), 8);

        graph.add_edge(n3, n3, EdgeData(500));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());

        graph.add_edge(n3, n3, EdgeData(500));
        assert!(graph.verify_node_pairing());
        assert!(!graph.verify_edge_mirror_property());
        assert_eq!(graph.edge_count(), 10);
    }
}
