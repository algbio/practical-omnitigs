use traitgraph::implementation::bit_vector_subgraph::BitVectorSubgraph;
use crate::restricted_reachability::{compute_hydrostructure_forward_reachability, compute_hydrostructure_backward_reachability};
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::interface::subgraph::Subgraph;
use traitgraph::walks::VecEdgeWalk;

pub enum Hydrostructure<Graph: GraphBase, SubgraphType> {
    BridgeLike {
        r_plus: SubgraphType,
        r_minus: SubgraphType,
        azb: VecEdgeWalk<Graph>,
    },
    Avertible {
        azb: VecEdgeWalk<Graph>,
    },
}

impl<'a, Graph: StaticGraph> Hydrostructure<Graph, BitVectorSubgraph<'a, Graph>> {
    pub fn compute_with_bitvector_subgraph(graph: &'a Graph, azb: VecEdgeWalk<Graph>) -> Self {
        let r_plus = compute_hydrostructure_forward_reachability(graph, &azb);
        let r_minus = compute_hydrostructure_backward_reachability(graph, &azb);

        if let (Some(r_plus), Some(r_minus)) = (r_plus, r_minus) {
            Self::BridgeLike {r_plus, r_minus, azb}
        } else {
            Self::Avertible {azb}
        }
    }
}

impl<'a, Graph: StaticGraph, SubgraphType: Subgraph<'a, Graph>> Hydrostructure<Graph, SubgraphType> {
    pub fn new_bridge_like(r_plus: SubgraphType, r_minus: SubgraphType, azb: VecEdgeWalk<Graph>) -> Self {
        Self::BridgeLike {
            r_plus, r_minus, azb
        }
    }
    pub fn new_avertible(azb: VecEdgeWalk<Graph>) -> Self {
        Self::Avertible {
            azb
        }
    }

    pub fn compute(graph: &'a Graph, azb: VecEdgeWalk<Graph>) -> Self {
        let r_plus = compute_hydrostructure_forward_reachability(graph, &azb);
        let r_minus = compute_hydrostructure_backward_reachability(graph, &azb);

        if let (Some(r_plus), Some(r_minus)) = (r_plus, r_minus) {
            Self::BridgeLike {r_plus, r_minus, azb}
        } else {
            Self::Avertible {azb}
        }
    }

    pub fn is_edge_river(&self, edge: Graph::EdgeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {r_plus, r_minus, azb} => !r_plus.contains_edge(edge) && !r_minus.contains_edge(edge),
            Hydrostructure::Avertible {azb} => false
        }
    }

    pub fn is_edge_vapor(&self, edge: Graph::EdgeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {r_plus, r_minus, azb} => r_plus.contains_edge(edge) && r_minus.contains_edge(edge),
            Hydrostructure::Avertible {azb} => true
        }
    }

    pub fn is_edge_cloud(&self, edge: Graph::EdgeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {r_plus, r_minus, azb} => !r_plus.contains_edge(edge) && r_minus.contains_edge(edge),
            Hydrostructure::Avertible {azb} => false
        }
    }

    pub fn is_edge_sea(&self, edge: Graph::EdgeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {r_plus, r_minus, azb} => r_plus.contains_edge(edge) && !r_minus.contains_edge(edge),
            Hydrostructure::Avertible {azb} => false
        }
    }

    pub fn is_node_river(&self, node: Graph::NodeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {r_plus, r_minus, azb} => !r_plus.contains_node(node) && !r_minus.contains_node(node),
            Hydrostructure::Avertible {azb} => false
        }
    }

    pub fn is_node_vapor(&self, node: Graph::NodeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {r_plus, r_minus, azb} => r_plus.contains_node(node) && r_minus.contains_node(node),
            Hydrostructure::Avertible {azb} => true
        }
    }

    pub fn is_node_cloud(&self, node: Graph::NodeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {r_plus, r_minus, azb} => !r_plus.contains_node(node) && r_minus.contains_node(node),
            Hydrostructure::Avertible {azb} => false
        }
    }

    pub fn is_node_sea(&self, node: Graph::NodeIndex) -> bool {
        match self {
            Hydrostructure::BridgeLike {r_plus, r_minus, azb} => r_plus.contains_node(node) && !r_minus.contains_node(node),
            Hydrostructure::Avertible {azb} => false
        }
    }
}