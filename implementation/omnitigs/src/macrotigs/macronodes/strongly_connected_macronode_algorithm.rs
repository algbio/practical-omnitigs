use super::{MacronodeAlgorithm, Macronodes};
use crate::unitigs::Unitigs;
use traitgraph::interface::StaticGraph;
use traitgraph::walks::EdgeWalk;

/// Compute the macronodes of a strongly connected graph.
pub struct StronglyConnectedMacronodes;

impl<Graph: StaticGraph> MacronodeAlgorithm<Graph> for StronglyConnectedMacronodes {
    fn compute_macronodes(graph: &Graph) -> Macronodes<Graph> {
        let unitigs = Unitigs::new(graph);
        let macronodes: Vec<_> = unitigs
            .into_iter()
            .filter(|unitig| {
                graph.out_degree(
                    graph
                        .edge_endpoints(unitig.iter().next().unwrap())
                        .from_node,
                ) == 1
                    && graph.in_degree(graph.edge_endpoints(unitig.iter().last().unwrap()).to_node)
                        == 1
            })
            .collect();
        Macronodes::new(macronodes)
    }
}
