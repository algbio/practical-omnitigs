use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecNodeWalk;

/// A macronode algorithm that requires the graph to be strongly connected.
pub mod strongly_connected_macronode_algorithm;

/// A struct containing the macronodes of an uncompressed graph, represented as walks through their uncompressed centers.
/// In an uncompressed graph, macronode centers are maximal unitigs with the property that their first node has outdegree = 1, and their last node has indegree = 1.
pub struct Macronodes<Graph: GraphBase> {
    macronodes: Vec<VecNodeWalk<Graph>>,
}

impl<Graph: GraphBase> Macronodes<Graph> {
    /// Creates a new `Macronodes` struct with the given vector of macronodes.
    pub fn new(macronodes: Vec<VecNodeWalk<Graph>>) -> Self {
        Self { macronodes }
    }

    /// Returns an iterator over the macronodes in this struct.
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = &'a VecNodeWalk<Graph>> {
        self.macronodes.iter()
    }
}

impl<'a, Graph: GraphBase> IntoIterator for &'a Macronodes<Graph> {
    type Item = &'a VecNodeWalk<Graph>;
    type IntoIter = std::slice::Iter<'a, VecNodeWalk<Graph>>;

    fn into_iter(self) -> Self::IntoIter {
        self.macronodes.iter()
    }
}

/// A trait abstracting over the concrete algorithm used to compute macronodes.
pub trait MacronodeAlgorithm<Graph: StaticGraph> {
    /// Compute all macronodes in a graph.
    fn compute_macronodes(graph: &Graph) -> Macronodes<Graph>;
}
