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
}
