use crate::hydrostructure::Hydrostructure;
use crate::restricted_reachability::{
    compute_hydrostructure_backward_reachability, compute_hydrostructure_forward_reachability,
};
use traitgraph::implementation::bit_vector_subgraph::BitVectorSubgraph;
use traitgraph::interface::subgraph::DecoratingSubgraph;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecEdgeWalk;

/// The hydrostructure for a walk aZb.
pub enum StaticHydrostructure<Graph: GraphBase, SubgraphType> {
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

impl<'a, Graph: StaticGraph> StaticHydrostructure<Graph, BitVectorSubgraph<'a, Graph>> {
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

impl<
        'a,
        Graph: 'a + StaticGraph,
        SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
    > StaticHydrostructure<Graph, SubgraphType>
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
    pub fn compute(graph: SubgraphType::ParentGraphRef, azb: VecEdgeWalk<Graph>) -> Self {
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

impl<
        'a,
        Graph: 'a + GraphBase,
        SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
    > Hydrostructure<Graph::NodeIndex, Graph::EdgeIndex>
    for StaticHydrostructure<Graph, SubgraphType>
{
    fn is_node_r_plus(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        match self {
            StaticHydrostructure::BridgeLike {
                r_plus,
                r_minus: _,
                azb: _,
            } => r_plus.contains_node(node),
            StaticHydrostructure::Avertible { azb: _ } => true,
        }
    }

    fn is_node_r_minus(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        match self {
            StaticHydrostructure::BridgeLike {
                r_plus: _,
                r_minus,
                azb: _,
            } => r_minus.contains_node(node),
            StaticHydrostructure::Avertible { azb: _ } => true,
        }
    }

    fn is_edge_r_plus(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        match self {
            StaticHydrostructure::BridgeLike {
                r_plus,
                r_minus: _,
                azb: _,
            } => r_plus.contains_edge(edge),
            StaticHydrostructure::Avertible { azb: _ } => true,
        }
    }

    fn is_edge_r_minus(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        match self {
            StaticHydrostructure::BridgeLike {
                r_plus: _,
                r_minus,
                azb: _,
            } => r_minus.contains_edge(edge),
            StaticHydrostructure::Avertible { azb: _ } => true,
        }
    }

    fn is_bridge_like(&self) -> bool {
        match self {
            StaticHydrostructure::BridgeLike {
                r_plus: _,
                r_minus: _,
                azb: _,
            } => true,
            StaticHydrostructure::Avertible { azb: _ } => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::StaticHydrostructure;
    use crate::hydrostructure::Hydrostructure;
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::subgraph::DecoratingSubgraph;
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
        let hydrostructure = StaticHydrostructure::compute_with_bitvector_subgraph(
            &graph,
            graph.create_edge_walk(&[e1, e2, e4]),
        );
        debug_assert!(hydrostructure.is_avertible());
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
        let hydrostructure = StaticHydrostructure::compute_with_bitvector_subgraph(
            &graph,
            graph.create_edge_walk(&[e1, e2, e3]),
        );
        debug_assert!(hydrostructure.is_avertible());
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
        let hydrostructure = StaticHydrostructure::compute_with_bitvector_subgraph(
            &graph,
            graph.create_edge_walk(&[e1, e2, e3]),
        );
        debug_assert!(hydrostructure.is_bridge_like());

        debug_assert!(hydrostructure.is_edge_cloud(e3));
        debug_assert!(hydrostructure.is_edge_sea(e1));
        debug_assert!(hydrostructure.is_edge_vapor(e2));
        debug_assert!(hydrostructure.is_edge_river(e5));

        match hydrostructure {
            StaticHydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => {
                debug_assert!(!r_plus.contains_node(n0));
                debug_assert!(r_plus.contains_node(n1));
                debug_assert!(r_plus.contains_node(n2));
                debug_assert!(!r_plus.contains_node(n3));
                debug_assert!(!r_plus.contains_node(n4));
                debug_assert!(!r_plus.contains_node(n5));

                debug_assert!(r_plus.contains_edge(e1));
                debug_assert!(r_plus.contains_edge(e2));
                debug_assert!(!r_plus.contains_edge(e3));
                debug_assert!(!r_plus.contains_edge(e4));
                debug_assert!(!r_plus.contains_edge(e5));
                debug_assert!(!r_plus.contains_edge(e6));
                debug_assert!(!r_plus.contains_edge(e7));
                debug_assert!(!r_plus.contains_edge(e8));

                debug_assert!(!r_minus.contains_node(n0));
                debug_assert!(r_minus.contains_node(n1));
                debug_assert!(r_minus.contains_node(n2));
                debug_assert!(!r_minus.contains_node(n3));
                debug_assert!(!r_minus.contains_node(n4));
                debug_assert!(!r_minus.contains_node(n5));

                debug_assert!(!r_minus.contains_edge(e1));
                debug_assert!(r_minus.contains_edge(e2));
                debug_assert!(r_minus.contains_edge(e3));
                debug_assert!(!r_minus.contains_edge(e4));
                debug_assert!(!r_minus.contains_edge(e5));
                debug_assert!(!r_minus.contains_edge(e6));
                debug_assert!(!r_minus.contains_edge(e7));
                debug_assert!(!r_minus.contains_edge(e8));
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
        let n6 = graph.add_node(6);
        let e1 = graph.add_edge(n0, n1, -1);
        let e2 = graph.add_edge(n1, n2, -2);
        let e3 = graph.add_edge(n2, n3, -3);
        let e4 = graph.add_edge(n3, n4, -4);
        let e5 = graph.add_edge(n3, n5, -5);
        let e6 = graph.add_edge(n4, n0, -6);
        let e7 = graph.add_edge(n5, n6, -7);
        let e8 = graph.add_edge(n1, n4, -8);
        let e9 = graph.add_edge(n5, n1, -9);
        let e10 = graph.add_edge(n6, n0, -10);
        let hydrostructure = StaticHydrostructure::compute_with_bitvector_subgraph(
            &graph,
            graph.create_edge_walk(&[e1, e2, e3]),
        );
        debug_assert!(hydrostructure.is_bridge_like());

        debug_assert!(hydrostructure.is_node_cloud(n3));
        debug_assert!(hydrostructure.is_node_sea(n0));
        debug_assert!(hydrostructure.is_node_vapor(n2));
        debug_assert!(hydrostructure.is_node_river(n6));
        debug_assert!(hydrostructure.is_edge_cloud(e3));
        debug_assert!(hydrostructure.is_edge_sea(e1));
        debug_assert!(hydrostructure.is_edge_vapor(e2));
        debug_assert!(hydrostructure.is_edge_river(e7));

        match hydrostructure {
            StaticHydrostructure::BridgeLike {
                r_plus,
                r_minus,
                azb: _,
            } => {
                debug_assert!(r_plus.contains_node(n0));
                debug_assert!(r_plus.contains_node(n1));
                debug_assert!(r_plus.contains_node(n2));
                debug_assert!(!r_plus.contains_node(n3));
                debug_assert!(r_plus.contains_node(n4));
                debug_assert!(!r_plus.contains_node(n5));
                debug_assert!(!r_plus.contains_node(n6));

                debug_assert!(r_plus.contains_edge(e1));
                debug_assert!(r_plus.contains_edge(e2));
                debug_assert!(!r_plus.contains_edge(e3));
                debug_assert!(!r_plus.contains_edge(e4));
                debug_assert!(!r_plus.contains_edge(e5));
                debug_assert!(r_plus.contains_edge(e6));
                debug_assert!(!r_plus.contains_edge(e7));
                debug_assert!(r_plus.contains_edge(e8));
                debug_assert!(!r_plus.contains_edge(e9));
                debug_assert!(!r_plus.contains_edge(e10));

                debug_assert!(!r_minus.contains_node(n0));
                debug_assert!(r_minus.contains_node(n1));
                debug_assert!(r_minus.contains_node(n2));
                debug_assert!(r_minus.contains_node(n3));
                debug_assert!(!r_minus.contains_node(n4));
                debug_assert!(r_minus.contains_node(n5));
                debug_assert!(!r_minus.contains_node(n6));

                debug_assert!(!r_minus.contains_edge(e1));
                debug_assert!(r_minus.contains_edge(e2));
                debug_assert!(r_minus.contains_edge(e3));
                debug_assert!(!r_minus.contains_edge(e4));
                debug_assert!(r_minus.contains_edge(e5));
                debug_assert!(!r_minus.contains_edge(e6));
                debug_assert!(!r_minus.contains_edge(e7));
                debug_assert!(!r_minus.contains_edge(e8));
                debug_assert!(r_minus.contains_edge(e9));
                debug_assert!(!r_minus.contains_edge(e10));
            }
            _ => panic!("Not bridge like"),
        }
    }
}
