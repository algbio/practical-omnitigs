use crate::restricted_reachability::{
    compute_hydrostructure_backward_reachability, compute_hydrostructure_forward_reachability,
};
use traitgraph::implementation::bit_vector_subgraph::BitVectorSubgraph;
use traitgraph::interface::subgraph::Subgraph;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecEdgeWalk;

/// The hydrostructure for a walk aZb.
pub enum Hydrostructure<Graph: GraphBase, SubgraphType> {
    /// In case the walk is bridge-like, `r_plus` and `r_minus` will be proper subgraphs.
    BridgeLike {
        /// The set `R⁺(aZb)`, defined as everything reachable from the first edge of `aZb` without using `aZb` as subwalk.
        r_plus: SubgraphType,
        /// The set `R⁻(aZb)`, defined as everything backwards reachable from the last edge of `aZb` without using `aZb` as subwalk.
        r_minus: SubgraphType,
        /// The walk the hydrostructure corresponds to.
        azb: VecEdgeWalk<Graph>,
    },
    /// In case the walk is avertible, the whole graph is in the vapor, so no subgraphs need to be stored.
    Avertible {
        /// The walk the hydrostructure corresponds to.
        azb: VecEdgeWalk<Graph>,
    },
}

impl<'a, Graph: StaticGraph> Hydrostructure<Graph, BitVectorSubgraph<'a, Graph>> {
    /// Compute the hydrostructure of a walk, representing `R⁺(aZb)` and `R⁻(aZb)` as `BitVectorSubgraph`s.
    pub fn compute_with_bitvector_subgraph(graph: &'a Graph, azb: VecEdgeWalk<Graph>) -> Self {
        let r_plus = compute_hydrostructure_forward_reachability(graph, &azb);
        let r_minus = compute_hydrostructure_backward_reachability(graph, &azb);

        if let (Some(r_plus), Some(r_minus)) = (r_plus, r_minus) {
            Self::BridgeLike {
                r_plus,
                r_minus,
                azb,
            }
        } else {
            Self::Avertible { azb }
        }
    }
}

impl<'a, Graph: StaticGraph, SubgraphType: Subgraph<'a, Graph>>
    Hydrostructure<Graph, SubgraphType>
{
    /// Initialise the hydrostructure of a _bridge-like_ walk `aZb` with given sets `R⁺(aZb)`, `R⁻(aZb)`.
    pub fn new_bridge_like(
        r_plus: SubgraphType,
        r_minus: SubgraphType,
        azb: VecEdgeWalk<Graph>,
    ) -> Self {
        Self::BridgeLike {
            r_plus,
            r_minus,
            azb,
        }
    }

    /// Initialise the hydrostructure of an _avertible_ walk `aZb`.
    pub fn new_avertible(azb: VecEdgeWalk<Graph>) -> Self {
        Self::Avertible { azb }
    }

    /// Compute the hydrostructure of a walk.
    pub fn compute(graph: &'a Graph, azb: VecEdgeWalk<Graph>) -> Self {
        let r_plus = compute_hydrostructure_forward_reachability(graph, &azb);
        let r_minus = compute_hydrostructure_backward_reachability(graph, &azb);

        if let (Some(r_plus), Some(r_minus)) = (r_plus, r_minus) {
            Self::BridgeLike {
                r_plus,
                r_minus,
                azb,
            }
        } else {
            Self::Avertible { azb }
        }
    }

    /// Returns true if the given edge is in the river.
    pub fn is_edge_river(&self, edge: Graph::EdgeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => !r_plus.contains_edge(edge) && !r_minus.contains_edge(edge),
            Hydrostructure::Avertible { azb: _ } => false,
        }
    }

    /// Returns true if the given edge is in the vapor.
    pub fn is_edge_vapor(&self, edge: Graph::EdgeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => r_plus.contains_edge(edge) && r_minus.contains_edge(edge),
            Hydrostructure::Avertible { azb: _ } => true,
        }
    }

    /// Returns true if the given edge is in the cloud.
    pub fn is_edge_cloud(&self, edge: Graph::EdgeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => !r_plus.contains_edge(edge) && r_minus.contains_edge(edge),
            Hydrostructure::Avertible { azb: _ } => false,
        }
    }

    /// Returns true if the given edge is in the sea.
    pub fn is_edge_sea(&self, edge: Graph::EdgeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => r_plus.contains_edge(edge) && !r_minus.contains_edge(edge),
            Hydrostructure::Avertible { azb: _ } => false,
        }
    }

    /// Returns true if the given node is in the river.
    pub fn is_node_river(&self, node: Graph::NodeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => !r_plus.contains_node(node) && !r_minus.contains_node(node),
            Hydrostructure::Avertible { azb: _ } => false,
        }
    }

    /// Returns true if the given node is in the vapor.
    pub fn is_node_vapor(&self, node: Graph::NodeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => r_plus.contains_node(node) && r_minus.contains_node(node),
            Hydrostructure::Avertible { azb: _ } => true,
        }
    }

    /// Returns true if the given node is in the cloud.
    pub fn is_node_cloud(&self, node: Graph::NodeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => !r_plus.contains_node(node) && r_minus.contains_node(node),
            Hydrostructure::Avertible { azb: _ } => false,
        }
    }

    /// Returns true if the given node is in the sea.
    pub fn is_node_sea(&self, node: Graph::NodeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => r_plus.contains_node(node) && !r_minus.contains_node(node),
            Hydrostructure::Avertible { azb: _ } => false,
        }
    }

    /// Returns true if the underlying walk of the hydrostructure is _bridge-like_.
    pub fn is_bridge_like(&self) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus: _,
                r_minus: _,
                azb: _,
            } => true,
            Hydrostructure::Avertible { azb: _ } => false,
        }
    }

    /// Returns true if the underlying walk of the hydrostructure is _avertible_.
    pub fn is_avertible(&self) -> bool {
        match self {
            Hydrostructure::BridgeLike {
                r_plus: _,
                r_minus: _,
                azb: _,
            } => false,
            Hydrostructure::Avertible { azb: _ } => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Hydrostructure;
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::subgraph::Subgraph;
    use traitgraph::interface::MutableGraphContainer;
    use traitgraph::interface::WalkableGraph;

    #[test]
    fn test_hydrostructure_avertible_by_shortcut() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let e1 = graph.add_edge(n0, n1, -1);
        let e2 = graph.add_edge(n1, n2, -2);
        graph.add_edge(n1, n2, -3);
        let e4 = graph.add_edge(n2, n3, -4);
        graph.add_edge(n3, n4, -5);
        graph.add_edge(n3, n5, -6);
        graph.add_edge(n4, n0, -7);
        graph.add_edge(n5, n0, -8);
        let hydrostructure = Hydrostructure::compute_with_bitvector_subgraph(
            &graph,
            graph.create_edge_walk(&[e1, e2, e4]),
        );
        assert!(hydrostructure.is_avertible());
    }

    #[test]
    fn test_hydrostructure_avertible_by_sea_cloud_edge() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let e1 = graph.add_edge(n0, n1, -1);
        let e2 = graph.add_edge(n1, n2, -2);
        let e3 = graph.add_edge(n2, n3, -3);
        graph.add_edge(n3, n4, -4);
        graph.add_edge(n3, n5, -5);
        graph.add_edge(n4, n0, -6);
        graph.add_edge(n5, n0, -7);
        graph.add_edge(n4, n2, -8);
        graph.add_edge(n1, n5, -9);
        graph.add_edge(n5, n4, -10);
        let hydrostructure = Hydrostructure::compute_with_bitvector_subgraph(
            &graph,
            graph.create_edge_walk(&[e1, e2, e3]),
        );
        assert!(hydrostructure.is_avertible());
    }

    #[test]
    fn test_hydrostructure_bridge_like_by_biunivocal() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let e1 = graph.add_edge(n0, n1, -1);
        let e2 = graph.add_edge(n1, n2, -2);
        let e3 = graph.add_edge(n2, n3, -3);
        let e4 = graph.add_edge(n3, n4, -4);
        let e5 = graph.add_edge(n3, n5, -5);
        let e6 = graph.add_edge(n4, n0, -6);
        let e7 = graph.add_edge(n5, n0, -7);
        let e8 = graph.add_edge(n5, n4, -8);
        let hydrostructure = Hydrostructure::compute_with_bitvector_subgraph(
            &graph,
            graph.create_edge_walk(&[e1, e2, e3]),
        );
        assert!(hydrostructure.is_bridge_like());
        match hydrostructure {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => {
                assert!(!r_plus.contains_node(n0));
                assert!(r_plus.contains_node(n1));
                assert!(r_plus.contains_node(n2));
                assert!(!r_plus.contains_node(n3));
                assert!(!r_plus.contains_node(n4));
                assert!(!r_plus.contains_node(n5));

                assert!(r_plus.contains_edge(e1));
                assert!(r_plus.contains_edge(e2));
                assert!(!r_plus.contains_edge(e3));
                assert!(!r_plus.contains_edge(e4));
                assert!(!r_plus.contains_edge(e5));
                assert!(!r_plus.contains_edge(e6));
                assert!(!r_plus.contains_edge(e7));
                assert!(!r_plus.contains_edge(e8));

                assert!(!r_minus.contains_node(n0));
                assert!(r_minus.contains_node(n1));
                assert!(r_minus.contains_node(n2));
                assert!(!r_minus.contains_node(n3));
                assert!(!r_minus.contains_node(n4));
                assert!(!r_minus.contains_node(n5));

                assert!(!r_minus.contains_edge(e1));
                assert!(r_minus.contains_edge(e2));
                assert!(r_minus.contains_edge(e3));
                assert!(!r_minus.contains_edge(e4));
                assert!(!r_minus.contains_edge(e5));
                assert!(!r_minus.contains_edge(e6));
                assert!(!r_minus.contains_edge(e7));
                assert!(!r_minus.contains_edge(e8));
            }
            _ => panic!("Not bridge like"),
        }
    }

    #[test]
    fn test_hydrostructure_bridge_like_non_trivial() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let e1 = graph.add_edge(n0, n1, -1);
        let e2 = graph.add_edge(n1, n2, -2);
        let e3 = graph.add_edge(n2, n3, -3);
        let e4 = graph.add_edge(n3, n4, -4);
        let e5 = graph.add_edge(n3, n5, -5);
        let e6 = graph.add_edge(n4, n0, -6);
        let e7 = graph.add_edge(n5, n0, -7);
        let e8 = graph.add_edge(n1, n4, -8);
        let e9 = graph.add_edge(n5, n1, -9);
        let hydrostructure = Hydrostructure::compute_with_bitvector_subgraph(
            &graph,
            graph.create_edge_walk(&[e1, e2, e3]),
        );
        assert!(hydrostructure.is_bridge_like());
        match hydrostructure {
            Hydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => {
                assert!(r_plus.contains_node(n0));
                assert!(r_plus.contains_node(n1));
                assert!(r_plus.contains_node(n2));
                assert!(!r_plus.contains_node(n3));
                assert!(r_plus.contains_node(n4));
                assert!(!r_plus.contains_node(n5));

                assert!(r_plus.contains_edge(e1));
                assert!(r_plus.contains_edge(e2));
                assert!(!r_plus.contains_edge(e3));
                assert!(!r_plus.contains_edge(e4));
                assert!(!r_plus.contains_edge(e5));
                assert!(r_plus.contains_edge(e6));
                assert!(!r_plus.contains_edge(e7));
                assert!(r_plus.contains_edge(e8));
                assert!(!r_plus.contains_edge(e9));

                assert!(!r_minus.contains_node(n0));
                assert!(r_minus.contains_node(n1));
                assert!(r_minus.contains_node(n2));
                assert!(r_minus.contains_node(n3));
                assert!(!r_minus.contains_node(n4));
                assert!(r_minus.contains_node(n5));

                assert!(!r_minus.contains_edge(e1));
                assert!(r_minus.contains_edge(e2));
                assert!(r_minus.contains_edge(e3));
                assert!(!r_minus.contains_edge(e4));
                assert!(r_minus.contains_edge(e5));
                assert!(!r_minus.contains_edge(e6));
                assert!(!r_minus.contains_edge(e7));
                assert!(!r_minus.contains_edge(e8));
                assert!(r_minus.contains_edge(e9));
            }
            _ => panic!("Not bridge like"),
        }
    }
}
