use super::{MacronodeAlgorithm, Macronodes};
use traitgraph::interface::StaticGraph;

/// Compute the macronodes of a strongly connected graph.
pub struct StronglyConnectedMacronodes;

impl<Graph: StaticGraph> MacronodeAlgorithm<Graph> for StronglyConnectedMacronodes {
    fn compute_macronodes(graph: &Graph) -> Macronodes<Graph> {
        todo!()
    }
}
