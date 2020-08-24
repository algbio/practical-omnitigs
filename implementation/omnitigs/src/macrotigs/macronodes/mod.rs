use crate::walks::VecWalk;
use traitgraph::interface::{GraphBase, StaticGraph};

/// A macronode algorithm that requires the graph to be strongly connected.
pub mod strongly_connected_macronode_algorithm;

/// A struct containing the macronodes of an uncompressed graph.
/// In an uncompressed graph, macronodes are maximal unitigs with the property that their first node has outdegree = 1, and their last node has indegree = 1.
pub struct Macronodes<Graph: GraphBase> {
    macronodes: Vec<VecWalk<Graph>>,
}

impl<Graph: GraphBase> Macronodes<Graph> {}

impl<'a, Graph: GraphBase> IntoIterator for &'a Macronodes<Graph> {
    type Item = &'a VecWalk<Graph>;
    type IntoIter = std::slice::Iter<'a, VecWalk<Graph>>;

    fn into_iter(self) -> Self::IntoIter {
        self.macronodes.iter()
    }
}

pub trait MacronodeAlgorithm<Graph: StaticGraph> {
    fn compute_macronodes(graph: &Graph) -> Macronodes<Graph>;
}
