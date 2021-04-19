use std::ops::{Index, Range};
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecNodeWalk;
use traitsequence::interface::Sequence;

/// A macronode algorithm that requires the graph to be strongly connected.
pub mod strongly_connected_macronode_algorithm;

/// A struct containing the macronodes of an uncompressed graph, represented as walks through their uncompressed centers.
/// In an uncompressed graph, macronode centers are maximal unitigs with the property that their first node has outdegree = 1, and their last node has indegree = 1.
pub struct Macronodes<Graph: GraphBase> {
    macronodes: Vec<VecNodeWalk<Graph>>,
}

impl<'a, Graph: GraphBase> Sequence<'a, VecNodeWalk<Graph>, [VecNodeWalk<Graph>]>
    for Macronodes<Graph>
where
    Graph::NodeIndex: 'a,
{
    type Iterator = std::slice::Iter<'a, VecNodeWalk<Graph>>;

    fn iter(&'a self) -> Self::Iterator {
        self.macronodes.iter()
    }

    fn len(&self) -> usize {
        self.macronodes.len()
    }
}

impl<Graph: GraphBase> Index<usize> for Macronodes<Graph> {
    type Output = VecNodeWalk<Graph>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.macronodes[index]
    }
}

impl<Graph: GraphBase> Index<Range<usize>> for Macronodes<Graph> {
    type Output = [VecNodeWalk<Graph>];

    fn index(&self, range: Range<usize>) -> &Self::Output {
        &self.macronodes[range]
    }
}

impl<Graph: GraphBase> From<Vec<VecNodeWalk<Graph>>> for Macronodes<Graph> {
    fn from(macronodes: Vec<VecNodeWalk<Graph>>) -> Self {
        Self { macronodes }
    }
}

/// A trait abstracting over the concrete algorithm used to compute macronodes.
pub trait MacronodeAlgorithm<Graph: StaticGraph> {
    /// Compute all macronodes in a graph.
    fn compute_macronodes(graph: &Graph) -> Macronodes<Graph>;
}
