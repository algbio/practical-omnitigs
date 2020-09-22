use crate::macrotigs::microtigs::Microtigs;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecEdgeWalk;

/// An algorithm to link maximal microtigs into maximal macrotigs.
pub mod default_macrotig_link_algorithm;

/// A structure containing macrotigs of a graph.
#[derive(Clone)]
pub struct Macrotigs<Graph: GraphBase> {
    macrotigs: Vec<VecEdgeWalk<Graph>>,
}

impl<Graph: GraphBase> Macrotigs<Graph> {
    /// Creates a new `Macrotigs` struct with the given vector of macrotigs.
    pub fn new(macrotigs: Vec<VecEdgeWalk<Graph>>) -> Self {
        Self { macrotigs }
    }

    /// Returns an iterator over the macrotigs in this struct.
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = &'a VecEdgeWalk<Graph>> {
        self.macrotigs.iter()
    }

    /// Returns the amount of macrotigs in this struct.
    pub fn len(&self) -> usize {
        self.macrotigs.len()
    }

    /// Returns true if this struct contains no macrotigs.
    pub fn is_empty(&self) -> bool {
        self.macrotigs.is_empty()
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for Macrotigs<Graph>
where
    Vec<VecEdgeWalk<Graph>>: std::ops::Index<IndexType>,
{
    type Output = <Vec<VecEdgeWalk<Graph>> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.macrotigs.index(index)
    }
}

impl<Graph: GraphBase> std::fmt::Debug for Macrotigs<Graph>
where
    Graph::NodeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Macrotigs[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase> PartialEq for Macrotigs<Graph>
where
    Graph::EdgeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.macrotigs == rhs.macrotigs
    }
}

impl<Graph: GraphBase> Eq for Macrotigs<Graph> where Graph::EdgeIndex: Eq {}

/// A trait abstracting over the concrete algorithm used to compute maximal macrotigs.
pub trait MaximalMacrotigsAlgorithm<Graph: StaticGraph> {
    /// Compute the maximal macrotigs from the given maximal microtigs.
    fn compute_maximal_macrotigs(graph: &Graph, microtigs: &Microtigs<Graph>) -> Macrotigs<Graph>;
}
