/// An algorithm to extract the maximal trivial omnitigs.
pub mod default_trivial_omnitigs;
/// An algorithm to extract non-trivial omnitigs from macrotigs using the incremental hydrostructure.
pub mod incremental_hydrostructure_macrotig_based_non_trivial_omnitigs;

use crate::macrotigs::macrotigs::Macrotigs;
use crate::omnitigs::default_trivial_omnitigs::DefaultTrivialOmnitigAlgorithm;
use crate::omnitigs::incremental_hydrostructure_macrotig_based_non_trivial_omnitigs::IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm;
use traitgraph::algo::traversal::univocal_traversal::univocal_extension_with_original_offset;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::{EdgeWalk, VecEdgeWalk};

/// An omnitig with information about its heart.
#[derive(Clone)]
pub struct Omnitig<Graph: GraphBase> {
    omnitig: VecEdgeWalk<Graph>,
    first_heart_edge: usize,
    last_heart_edge: usize,
}

impl<Graph: StaticGraph> Omnitig<Graph> {
    /// Computes the omnitig from the given omnitig heart.
    /// Does not check if the given walk is actually an omnitig heart.
    pub fn compute_from_heart(graph: &Graph, heart: &[Graph::EdgeIndex]) -> Self {
        let (first_heart_edge, univocal_extension) =
            univocal_extension_with_original_offset(graph, heart);
        let last_heart_edge = first_heart_edge + heart.len() - 1;
        Self::new(univocal_extension, first_heart_edge, last_heart_edge)
    }

    /// Computes the omnitig from the given superwalk of an non-trivial omnitig heart.
    /// The superwalk must still be a subwalk of the omnitig.
    /// Does not check if the given walk is actually an omnitig heart.
    /// Panics if the superwalk of the non-trivial omnitig heart does not have its first join edge before its last split edge.
    pub fn compute_from_non_trivial_heart_superwalk(
        graph: &Graph,
        heart_superwalk: &[Graph::EdgeIndex],
    ) -> Self {
        let mut omnitig = Self::compute_from_heart(graph, heart_superwalk);
        while !graph.is_join_edge(omnitig[omnitig.first_heart_edge]) {
            omnitig.first_heart_edge += 1;
            assert!(
                omnitig.first_heart_edge < omnitig.last_heart_edge,
                "First join is not before last split"
            );
        }
        while !graph.is_split_edge(omnitig[omnitig.last_heart_edge]) {
            omnitig.last_heart_edge -= 1;
            assert!(
                omnitig.first_heart_edge < omnitig.last_heart_edge,
                "First join is not before last split"
            );
        }
        omnitig
    }
}

impl<Graph: GraphBase> Omnitig<Graph> {
    /// Construct an `Omnitig` with the given attributes.
    pub fn new(edges: VecEdgeWalk<Graph>, first_heart_edge: usize, last_heart_edge: usize) -> Self {
        Self {
            omnitig: edges,
            first_heart_edge,
            last_heart_edge,
        }
    }

    /// Returns an iterator over the edges in the heart of this omnitig.
    pub fn iter_heart<'a>(&'a self) -> impl 'a + Iterator<Item = Graph::EdgeIndex> {
        self.omnitig
            .iter()
            .take(self.last_heart_edge + 1)
            .skip(self.first_heart_edge)
    }

    /// Returns a slice of the heart edges of this omnitig.
    pub fn heart(&self) -> &[Graph::EdgeIndex] {
        &self.omnitig[self.first_heart_edge..=self.last_heart_edge]
    }

    /// Returns the amount of omnitigs in this struct.
    pub fn len_heart(&self) -> usize {
        self.heart().len()
    }
}

impl<'a, Graph: GraphBase> EdgeWalk<'a, Graph> for Omnitig<Graph>
where
    Graph::EdgeIndex: 'a,
{
    type Iter = std::iter::Cloned<std::slice::Iter<'a, Graph::EdgeIndex>>;

    fn iter(&'a self) -> Self::Iter {
        self.omnitig.iter()
    }
}

impl<Graph: GraphBase> PartialEq for Omnitig<Graph>
where
    Graph::EdgeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.omnitig == rhs.omnitig
            && self.first_heart_edge == rhs.first_heart_edge
            && self.last_heart_edge == rhs.last_heart_edge
    }
}

impl<Graph: GraphBase> Eq for Omnitig<Graph> where Graph::EdgeIndex: Eq {}

impl<Graph: GraphBase> std::fmt::Debug for Omnitig<Graph>
where
    Graph::EdgeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Omnitig[")?;
        if let Some((i, first)) = self.iter().enumerate().next() {
            if i == self.first_heart_edge {
                write!(f, "|")?;
            }
            write!(f, "{:?}", first)?;
            if i == self.last_heart_edge {
                write!(f, "|")?;
            }
        }
        for (i, edge) in self.iter().enumerate().skip(1) {
            write!(f, ", ")?;
            if i == self.first_heart_edge {
                write!(f, "|")?;
            }
            write!(f, "{:?}", edge)?;
            if i == self.last_heart_edge {
                write!(f, "|")?;
            }
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for Omnitig<Graph>
where
    Vec<Graph::EdgeIndex>: std::ops::Index<IndexType>,
{
    type Output = <Vec<Graph::EdgeIndex> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.omnitig.index(index)
    }
}

/// A structure containing omnitigs of a graph.
#[derive(Clone, Default)]
pub struct Omnitigs<Graph: GraphBase> {
    omnitigs: Vec<Omnitig<Graph>>,
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
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = &'a Omnitig<Graph>> {
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
    pub fn push(&mut self, omnitig: Omnitig<Graph>) {
        self.omnitigs.push(omnitig);
    }
}

impl<Graph: GraphBase> From<Vec<Omnitig<Graph>>> for Omnitigs<Graph> {
    fn from(omnitigs: Vec<Omnitig<Graph>>) -> Self {
        Self { omnitigs }
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for Omnitigs<Graph>
where
    Vec<Omnitig<Graph>>: std::ops::Index<IndexType>,
{
    type Output = <Vec<Omnitig<Graph>> as std::ops::Index<IndexType>>::Output;

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
    use crate::omnitigs::{Omnitig, Omnitigs};
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
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e4, e20]), 2, 3),
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e5, e21]), 2, 3),
                Omnitig::new(graph.create_edge_walk(&[e28, e27, e26, e6, e0, e1]), 0, 3),
                Omnitig::new(graph.create_edge_walk(&[e14, e7, e0, e1]), 0, 1),
                Omnitig::new(graph.create_edge_walk(&[e15, e8, e0, e1]), 0, 1),
                Omnitig::new(graph.create_edge_walk(&[e16, e9, e0, e1]), 0, 1),
                Omnitig::new(
                    graph.create_edge_walk(&[e0, e1, e2, e10, e22, e23, e24]),
                    3,
                    6
                ),
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e2, e11, e17]), 3, 4),
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e3, e12, e18]), 3, 4),
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e3, e13, e19]), 3, 4),
                Omnitig::new(graph.create_edge_walk(&[e25]), 0, 0),
                Omnitig::new(graph.create_edge_walk(&[e29]), 0, 0),
            ])
        );
    }
}
