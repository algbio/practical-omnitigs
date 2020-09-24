/// An algorithm to extract the maximal trivial omnitigs.
pub mod default_trivial_omnitigs;
/// An algorithm to extract non-trivial omnitigs from macrotigs using the incremental hydrostructure.
pub mod incremental_hydrostructure_macrotig_based_non_trivial_omnitigs;

use crate::macrotigs::macrotigs::Macrotigs;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecEdgeWalk;

/// A structure containing omnitigs of a graph.
#[derive(Clone)]
pub struct Omnitigs<Graph: GraphBase> {
    omnitigs: Vec<VecEdgeWalk<Graph>>,
}

impl<Graph: GraphBase> Omnitigs<Graph> {
    /// Creates a new `Omnitigs` struct with the given vector of omnitigs.
    pub fn new(omnitigs: Vec<VecEdgeWalk<Graph>>) -> Self {
        Self { omnitigs }
    }

    /// Returns an iterator over the omnitigs in this struct.
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = &'a VecEdgeWalk<Graph>> {
        self.omnitigs.iter()
    }

    /// Returns the amount of omnitigs in this struct.
    pub fn len(&self) -> usize {
        self.omnitigs.len()
    }

    /// Returns true if this struct contains no omnitigs.
    pub fn is_empty(&self) -> bool {
        self.omnitigs.is_empty()
    }

    /// Adds the given omnitig to this struct.
    pub fn push(&mut self, omnitig: VecEdgeWalk<Graph>) {
        self.omnitigs.push(omnitig);
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for Omnitigs<Graph>
where
    Vec<VecEdgeWalk<Graph>>: std::ops::Index<IndexType>,
{
    type Output = <Vec<VecEdgeWalk<Graph>> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.omnitigs.index(index)
    }
}

impl<Graph: GraphBase> std::fmt::Debug for Omnitigs<Graph>
where
    Graph::NodeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Omnitigs[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase> PartialEq for Omnitigs<Graph>
where
    Graph::EdgeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.omnitigs == rhs.omnitigs
    }
}

impl<Graph: GraphBase> Eq for Omnitigs<Graph> where Graph::EdgeIndex: Eq {}

/// A trait abstracting over the concrete algorithm used to compute maximal non-trivial omnitigs based on macrotigs.
pub trait MacrotigBasedNonTrivialOmnitigAlgorithm<Graph: StaticGraph> {
    /// Compute the maximal non-trivial omnitigs of the given the maximal macrotigs.
    fn compute_maximal_non_trivial_omnitigs(
        graph: &Graph,
        macrotigs: &Macrotigs<Graph>,
    ) -> Omnitigs<Graph>;
}

/// A trait abstracting over the concrete algorithm used to compute maximal trivial omnitigs.
pub trait TrivialOmnitigAlgorithm<Graph: StaticGraph> {
    /// To a sequence of maximal non-trivial omnitigs add the maximal trivial omnitigs.
    /// The function should not compute any trivial omnitigs that are subwalks of maximal non-trivial omnitigs.
    fn compute_maximal_trivial_omnitigs(
        graph: &Graph,
        omnitigs: Omnitigs<Graph>,
    ) -> Omnitigs<Graph>;
}