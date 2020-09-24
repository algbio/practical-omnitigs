/// An algorithm to extract the maximal trivial omnitigs.
pub mod default_trivial_omnitigs;
/// An algorithm to extract non-trivial omnitigs from macrotigs using the incremental hydrostructure.
pub mod incremental_hydrostructure_macrotig_based_non_trivial_omnitigs;

use crate::macrotigs::macrotigs::Macrotigs;
use crate::omnitigs::default_trivial_omnitigs::DefaultTrivialOmnitigAlgorithm;
use crate::omnitigs::incremental_hydrostructure_macrotig_based_non_trivial_omnitigs::IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecEdgeWalk;

/// A structure containing omnitigs of a graph.
#[derive(Clone, Default)]
pub struct Omnitigs<Graph: GraphBase> {
    omnitigs: Vec<VecEdgeWalk<Graph>>,
}

impl<Graph: StaticGraph> Omnitigs<Graph> {
    /// Computes the maximal omnitigs of the given graph.
    pub fn compute(graph: &Graph) -> Self {
        let maximal_macrotigs = Macrotigs::compute(graph);
        let maximal_non_trivial_omnitigs = IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm::compute_maximal_non_trivial_omnitigs(graph, &maximal_macrotigs);
        DefaultTrivialOmnitigAlgorithm::compute_maximal_trivial_omnitigs(
            graph,
            maximal_non_trivial_omnitigs,
        )
    }

    /// Computes the maximal trivial omnitigs of the given graph, including those that are subwalks of maximal non-trivial omnitigs.
    pub fn compute_trivial_only(graph: &Graph) -> Self {
        DefaultTrivialOmnitigAlgorithm::compute_maximal_trivial_omnitigs(graph, Omnitigs::new())
    }
}

impl<Graph: GraphBase> Omnitigs<Graph> {
    /// Creates a new empty `Omnitigs` struct.
    pub fn new() -> Self {
        Self {
            omnitigs: Default::default(),
        }
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

impl<Graph: GraphBase> From<Vec<VecEdgeWalk<Graph>>> for Omnitigs<Graph> {
    fn from(omnitigs: Vec<VecEdgeWalk<Graph>>) -> Self {
        Self { omnitigs }
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

#[cfg(test)]
mod tests {
    use crate::omnitigs::Omnitigs;
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::MutableGraphContainer;
    use traitgraph::interface::WalkableGraph;

    #[test]
    fn test_compute_only_trivial_omnitigs_simple() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let n4 = graph.add_node(());
        let n5 = graph.add_node(());
        let n6 = graph.add_node(());
        let n7 = graph.add_node(());
        let n8 = graph.add_node(());
        let n9 = graph.add_node(());
        let n10 = graph.add_node(());
        let n11 = graph.add_node(());
        let n12 = graph.add_node(());
        let n13 = graph.add_node(());
        let n14 = graph.add_node(());
        let n15 = graph.add_node(());
        let n16 = graph.add_node(());
        let n17 = graph.add_node(());
        let n18 = graph.add_node(());
        let n19 = graph.add_node(());
        let n20 = graph.add_node(());

        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n2, n3, ());
        let e3 = graph.add_edge(n2, n4, ());
        let e4 = graph.add_edge(n2, n5, ());
        let e5 = graph.add_edge(n2, n6, ());
        let e6 = graph.add_edge(n7, n0, ()); // Comes from all except n11.
        let e7 = graph.add_edge(n8, n0, ());
        let e8 = graph.add_edge(n9, n0, ());
        let e9 = graph.add_edge(n10, n0, ());
        let e10 = graph.add_edge(n3, n11, ()); // Goes to all except n7.
        let e11 = graph.add_edge(n3, n12, ());
        let e12 = graph.add_edge(n4, n13, ());
        let e13 = graph.add_edge(n4, n14, ());
        let e14 = graph.add_edge(n17, n8, ());
        let e15 = graph.add_edge(n17, n9, ());
        let e16 = graph.add_edge(n17, n10, ());
        let e17 = graph.add_edge(n12, n18, ());
        let e18 = graph.add_edge(n13, n18, ());
        let e19 = graph.add_edge(n14, n18, ());
        let e20 = graph.add_edge(n5, n18, ());
        let e21 = graph.add_edge(n6, n18, ());
        let e22 = graph.add_edge(n11, n15, ());
        let e23 = graph.add_edge(n15, n16, ());
        let e24 = graph.add_edge(n16, n17, ());
        let e25 = graph.add_edge(n17, n17, ());
        let e26 = graph.add_edge(n20, n7, ());
        let e27 = graph.add_edge(n19, n20, ());
        let e28 = graph.add_edge(n18, n19, ());
        let e29 = graph.add_edge(n18, n18, ());

        let maximal_trivial_omnitigs = Omnitigs::compute_trivial_only(&graph);
        assert_eq!(
            maximal_trivial_omnitigs,
            Omnitigs::from(vec![
                graph.create_edge_walk(&[e0, e1, e4, e20]),
                graph.create_edge_walk(&[e0, e1, e5, e21]),
                graph.create_edge_walk(&[e28, e27, e26, e6, e0, e1]),
                graph.create_edge_walk(&[e14, e7, e0, e1]),
                graph.create_edge_walk(&[e15, e8, e0, e1]),
                graph.create_edge_walk(&[e16, e9, e0, e1]),
                graph.create_edge_walk(&[e0, e1, e2, e10, e22, e23, e24]),
                graph.create_edge_walk(&[e0, e1, e2, e11, e17]),
                graph.create_edge_walk(&[e0, e1, e3, e12, e18]),
                graph.create_edge_walk(&[e0, e1, e3, e13, e19]),
                graph.create_edge_walk(&[e25]),
                graph.create_edge_walk(&[e29]),
            ])
        );
    }
}
