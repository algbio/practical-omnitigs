use crate::macrotigs::macronodes::Macronodes;
use super::{MaximalMicrotigsAlgorithm, Microtigs};
use traitgraph::interface::StaticGraph;
use traitgraph::walks::VecEdgeWalk;
use traitgraph::walks::NodeWalk;
use traitgraph::implementation::bit_vector_subgraph::BitVectorSubgraph;
use crate::restricted_reachability::compute_restricted_backward_reachability;

/// Compute the maximal microtigs of a strongly connected graph using hydrostructure-based queries.
pub struct StronglyConnectedHydrostructureBasedMaximalMicrotigs;

impl<Graph: StaticGraph> MaximalMicrotigsAlgorithm<Graph> for StronglyConnectedHydrostructureBasedMaximalMicrotigs
    where Graph::EdgeIndex: 'static
{
    fn compute_maximal_microtigs(graph: &Graph, macronodes: &Macronodes<Graph>) -> Microtigs<Graph> {
        let mut result = Vec::new();
        todo!();

        for macronode in macronodes {
            assert!(!macronode.is_empty());

            //////////////////////////////////////////
            //// Compute right-maximal microtigs. ////
            //////////////////////////////////////////
            let rmm1 = Vec::new();
            // Choose any incoming edge of the first node of the macronode center, and compute its
            // restricted backwards reachability.
            let first_edge = graph.in_neighbors(macronode.iter().next().unwrap()).next().expect("Found sink, but this algorithm requires the graph to be strongly connected").edge_id;
            rmm1.push(first_edge);
            let r_minus: BitVectorSubgraph<_> = compute_restricted_backward_reachability(graph, first_edge);

            // Compute left-maximal microtigs.

            // Stitch them together.
            // TODO Here one could also stitch the two maximal microtigs together if they meet in a self-bivalent edge. But we do not do that for now.
        }

        Microtigs::new(result)
    }
}