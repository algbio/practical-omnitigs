use crate::macrotigs::microtigs::{Microtigs, MaximalMicrotigsAlgorithm};
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecEdgeWalk;
use crate::macrotigs::macronodes::strongly_connected_macronode_algorithm::StronglyConnectedMacronodes;
use crate::macrotigs::microtigs::strongly_connected_hydrostructure_based_maximal_microtig_algorithm::StronglyConnectedHydrostructureBasedMaximalMicrotigs;
use crate::macrotigs::macronodes::MacronodeAlgorithm;
use crate::macrotigs::macrotigs::default_macrotig_link_algorithm::DefaultMacrotigLinkAlgorithm;
use traitsequence::interface::Sequence;

/// An algorithm to link maximal microtigs into maximal macrotigs.
pub mod default_macrotig_link_algorithm;

/// A structure containing macrotigs of a graph.
#[derive(Clone, Default)]
pub struct Macrotigs<Graph: GraphBase> {
    macrotigs: Vec<VecEdgeWalk<Graph>>,
}

impl<Graph: StaticGraph> Macrotigs<Graph> {
    /// Computes the maximal macrotigs of the given graph.
    pub fn compute(graph: &Graph) -> Self {
        let macronodes = StronglyConnectedMacronodes::compute_macronodes(graph);
        let maximal_microtigs =
            StronglyConnectedHydrostructureBasedMaximalMicrotigs::compute_maximal_microtigs(
                graph,
                &macronodes,
            );
        DefaultMacrotigLinkAlgorithm::compute_maximal_macrotigs(graph, &maximal_microtigs)
    }
}

impl<Graph: GraphBase> Macrotigs<Graph> {
    /// Creates a new empty `Macrotigs` struct.
    pub fn new() -> Self {
        Self {
            macrotigs: Default::default(),
        }
    }
}

impl<'a, Graph: 'a + GraphBase> Sequence<'a, VecEdgeWalk<Graph>, [VecEdgeWalk<Graph>]>
    for Macrotigs<Graph>
{
    type Iterator = std::slice::Iter<'a, VecEdgeWalk<Graph>>;

    fn iter(&'a self) -> Self::Iterator {
        self.macrotigs.iter()
    }

    fn len(&self) -> usize {
        self.macrotigs.len()
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

impl<Graph: GraphBase> From<Vec<VecEdgeWalk<Graph>>> for Macrotigs<Graph> {
    fn from(macrotigs: Vec<VecEdgeWalk<Graph>>) -> Self {
        Self { macrotigs }
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
