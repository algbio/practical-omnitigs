use crate::macrotigs::macronodes::Macronodes;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecEdgeWalk;

/// A maximal microtig algorithm that uses hydrostructure-based queries and requires the graph to be strongly connected.
pub mod strongly_connected_hydrostructure_based_maximal_microtig_algorithm;

/// A structure containing microtigs of a graph.
///
/// Since we do not compress our graph, the bivalent edges that cause a microtig to end might be paths.
/// To make connecting microtigs simpler, we therefore define a maximal microtig as ending and starting
/// with the last edge of a _bivalent path_.
#[derive(Clone)]
pub struct Microtigs<Graph: GraphBase> {
    microtigs: Vec<VecEdgeWalk<Graph>>,
}

impl<Graph: GraphBase> Microtigs<Graph> {
    /// Creates a new `Microtigs` struct with the given vector of microtigs.
    pub fn new(microtigs: Vec<VecEdgeWalk<Graph>>) -> Self {
        Self { microtigs }
    }

    /// Returns an iterator over the microtigs in this struct.
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = &'a VecEdgeWalk<Graph>> {
        self.microtigs.iter()
    }

    /// Returns the amount of microtigs in this struct.
    pub fn len(&self) -> usize {
        self.microtigs.len()
    }

    /// Returns true if this struct contains no microtigs.
    pub fn is_empty(&self) -> bool {
        self.microtigs.is_empty()
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for Microtigs<Graph>
where
    Vec<VecEdgeWalk<Graph>>: std::ops::Index<IndexType>,
{
    type Output = <Vec<VecEdgeWalk<Graph>> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.microtigs.index(index)
    }
}

impl<Graph: GraphBase> std::fmt::Debug for Microtigs<Graph>
where
    Graph::NodeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Microtigs[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase> PartialEq for Microtigs<Graph>
where
    Graph::EdgeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.microtigs == rhs.microtigs
    }
}

impl<Graph: GraphBase> Eq for Microtigs<Graph> where Graph::EdgeIndex: Eq {}

/// A trait abstracting over the concrete algorithm used to compute maximal microtigs.
pub trait MaximalMicrotigsAlgorithm<Graph: StaticGraph> {
    /// Compute the maximal microtigs of the given macronode centers.
    fn compute_maximal_microtigs(graph: &Graph, macronodes: &Macronodes<Graph>)
        -> Microtigs<Graph>;
}
