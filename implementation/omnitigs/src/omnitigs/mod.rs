/// An algorithm to extract the maximal node-centric trivial omnitigs.
pub mod default_node_centric_trivial_omnitigs;
/// An algorithm to extract the maximal trivial omnitigs.
pub mod default_trivial_omnitigs;
/// An algorithm to extract non-trivial omnitigs from macrotigs using the incremental hydrostructure.
pub mod incremental_hydrostructure_macrotig_based_non_trivial_omnitigs;
/// Different algorithms to compute univocal extensions.
pub mod univocal_extension_algorithms;

use crate::macrotigs::macrotigs::Macrotigs;
use crate::omnitigs::default_node_centric_trivial_omnitigs::DefaultTrivialNodeCentricOmnitigAlgorithm;
use crate::omnitigs::default_trivial_omnitigs::{
    NonSccTrivialOmnitigAlgorithm, SccTrivialOmnitigAlgorithm,
};
use crate::omnitigs::incremental_hydrostructure_macrotig_based_non_trivial_omnitigs::IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm;
use crate::omnitigs::univocal_extension_algorithms::{
    NonSccNodeCentricUnivocalExtensionStrategy, SccNodeCentricUnivocalExtensionStrategy,
};
use crate::walks::EdgeOmnitigLikeExt;
use bigraph::interface::static_bigraph::{StaticBigraph, StaticEdgeCentricBigraph};
use bigraph::interface::BidirectedData;
use std::cmp::Ordering;
use std::iter::FromIterator;
use traitgraph::index::GraphIndex;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::{EdgeWalk, VecEdgeWalk, VecNodeWalk};
use traitsequence::interface::Sequence;

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
            heart.compute_univocal_extension_with_original_offset(graph);
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
            debug_assert!(
                omnitig.first_heart_edge < omnitig.last_heart_edge,
                "First join is not before last split"
            );
        }
        while !graph.is_split_edge(omnitig[omnitig.last_heart_edge]) {
            omnitig.last_heart_edge -= 1;
            debug_assert!(
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
    pub fn iter_heart(&self) -> impl Iterator<Item = &Graph::EdgeIndex> {
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

impl<'a, Graph: GraphBase> EdgeWalk<'a, Graph, [Graph::EdgeIndex]> for Omnitig<Graph> where
    Graph::EdgeIndex: 'a
{
}

impl<'a, Graph: GraphBase> Sequence<'a, Graph::EdgeIndex, [Graph::EdgeIndex]> for Omnitig<Graph>
where
    Graph::EdgeIndex: 'a,
{
    type Iterator =
        <VecEdgeWalk<Graph> as Sequence<'a, Graph::EdgeIndex, [Graph::EdgeIndex]>>::Iterator;

    fn iter(&'a self) -> Self::Iterator {
        self.omnitig.iter()
    }

    fn len(&self) -> usize {
        self.omnitig.len()
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
#[derive(Clone)]
pub struct Omnitigs<Graph: GraphBase> {
    omnitigs: Vec<Omnitig<Graph>>,
    omnitigs_per_macrotig: Vec<usize>,
}

impl<Graph: StaticGraph> Omnitigs<Graph> {
    /// Computes the maximal omnitigs of the given graph.
    pub fn compute(graph: &Graph) -> Self {
        let maximal_macrotigs = Macrotigs::compute(graph);
        let maximal_non_trivial_omnitigs = IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm::compute_maximal_non_trivial_omnitigs(graph, &maximal_macrotigs);
        SccTrivialOmnitigAlgorithm::compute_maximal_trivial_omnitigs(
            graph,
            maximal_non_trivial_omnitigs,
        )
    }

    /// Computes the maximal non-trivial omnitigs of the graph.
    pub fn compute_non_trivial_only(graph: &Graph) -> Self {
        let maximal_macrotigs = Macrotigs::compute(graph);
        IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm::compute_maximal_non_trivial_omnitigs(graph, &maximal_macrotigs)
    }

    /// Computes the maximal trivial omnitigs of the given graph, including those that are subwalks of maximal non-trivial omnitigs.
    pub fn compute_trivial_only(graph: &Graph) -> Self {
        SccTrivialOmnitigAlgorithm::compute_maximal_trivial_omnitigs(graph, Omnitigs::default())
    }

    /// Computes the maximal trivial omnitigs of the given graph, including those that are subwalks of maximal non-trivial omnitigs.
    /// This algorithm allows the graph to be not strongly connected, but it is a bit slower, especially for long trivial omnitigs.
    pub fn compute_trivial_only_non_scc(graph: &Graph) -> Self {
        NonSccTrivialOmnitigAlgorithm::compute_maximal_trivial_omnitigs(graph, Omnitigs::default())
    }
}

impl<Graph: StaticEdgeCentricBigraph> Omnitigs<Graph>
where
    Graph::EdgeData: BidirectedData + Eq,
    Graph::NodeData: std::fmt::Debug,
{
    /// Retains only one direction of each pair of reverse-complemental omnitigs.
    ///
    /// Note: I am not sure if this method is correct in all cases, but it will panic if it finds a case where it is not correct.
    ///       For practical genomes it seems to work.
    pub fn remove_reverse_complements(&mut self, graph: &Graph) {
        // Maps from edges to omnitigs that have this edge as first edge in their heart.
        let mut first_heart_edge_map = vec![usize::max_value(); graph.edge_count()];
        for (i, omnitig) in self.iter().enumerate() {
            let first_heart_edge = omnitig.iter_heart().next().expect("Omnitig has no heart");
            // I am not sure if the following assumption is correct.
            debug_assert_eq!(
                first_heart_edge_map[first_heart_edge.as_usize()],
                usize::max_value(),
                "Found two omnitigs hearts starting with the same edge."
            );
            first_heart_edge_map[first_heart_edge.as_usize()] = i;
        }

        let mut retain_indices = Vec::with_capacity(self.len());
        for (i, omnitig) in self.iter().enumerate() {
            let reverse_complement_first_heart_edge = graph
                .mirror_edge_edge_centric(
                    *omnitig.iter_heart().last().expect("Omnitig has no heart."),
                )
                .expect("Edge has no reverse complement.");
            let reverse_complement_candidate_index =
                first_heart_edge_map[reverse_complement_first_heart_edge.as_usize()];
            if reverse_complement_candidate_index < i {
                let reverse_complement_candidate = &self[reverse_complement_candidate_index];
                for (edge, reverse_complement_edge) in omnitig
                    .iter()
                    .zip(reverse_complement_candidate.iter().rev())
                {
                    let complements_complement_edge = graph
                        .mirror_edge_edge_centric(*reverse_complement_edge)
                        .expect("Edge has no reverse complement.");
                    // If our algorithms are sound, then this assumption should be correct.
                    debug_assert_eq!(
                        *edge,
                        complements_complement_edge,
                        "Found reverse complement candidate, but it is not a reverse complement:\nomnitig: {:?}\nnode omnitig: {:?}\nomnitig indegree:  {}\nomnitig outdegree: {}\nrevcomp: {:?}\nnode revcomp: {:?}\nrevcomp indegree:  {}\nrevcomp outdegree: {}",
                        omnitig,
                        omnitig.clone_as_node_walk::<VecNodeWalk<Graph>>(graph).unwrap().iter().map(|&n| graph.node_data(n)).collect::<Vec<_>>(),
                        graph.in_degree(*omnitig.clone_as_node_walk::<VecNodeWalk<Graph>>(graph).unwrap().first().unwrap()),
                        graph.out_degree(*omnitig.clone_as_node_walk::<VecNodeWalk<Graph>>(graph).unwrap().last().unwrap()),
                        reverse_complement_candidate,
                        reverse_complement_candidate.clone_as_node_walk::<VecNodeWalk<Graph>>(graph).unwrap().iter().map(|&n| graph.node_data(n)).collect::<Vec<_>>(),
                        graph.in_degree(*reverse_complement_candidate.clone_as_node_walk::<VecNodeWalk<Graph>>(graph).unwrap().first().unwrap()),
                        graph.out_degree(*reverse_complement_candidate.clone_as_node_walk::<VecNodeWalk<Graph>>(graph).unwrap().last().unwrap()),
                    );
                }
            } else {
                retain_indices.push(i);
            }
        }

        let mut omnitigs = Vec::new();
        std::mem::swap(&mut omnitigs, &mut self.omnitigs);
        for (i, omnitig) in omnitigs.into_iter().enumerate() {
            if self.omnitigs.len() == retain_indices.len() {
                break;
            }

            if i == retain_indices[self.omnitigs.len()] {
                self.omnitigs.push(omnitig);
            }
        }
    }
}

impl<Graph: GraphBase> Omnitigs<Graph> {
    /// Creates a new `Omnitigs` struct from the given omnitigs and statistics.
    pub fn new(omnitigs: Vec<Omnitig<Graph>>, omnitigs_per_macrotig: Vec<usize>) -> Self {
        Self {
            omnitigs,
            omnitigs_per_macrotig,
        }
    }

    /// Returns an iterator over the omnitigs in this struct.
    pub fn iter(&self) -> impl Iterator<Item = &Omnitig<Graph>> {
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

    /// Returns a slice of omnitig counts per macrotig.
    pub fn omnitigs_per_macrotig(&self) -> &[usize] {
        &self.omnitigs_per_macrotig
    }
}

impl<Graph: GraphBase> Default for Omnitigs<Graph> {
    fn default() -> Self {
        Self {
            omnitigs: Default::default(),
            omnitigs_per_macrotig: Default::default(),
        }
    }
}

impl<Graph: GraphBase> From<Vec<Omnitig<Graph>>> for Omnitigs<Graph> {
    fn from(omnitigs: Vec<Omnitig<Graph>>) -> Self {
        Self {
            omnitigs,
            omnitigs_per_macrotig: Default::default(),
        }
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

impl<Graph: GraphBase> IntoIterator for Omnitigs<Graph> {
    type Item = Omnitig<Graph>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.omnitigs.into_iter()
    }
}

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
    /// The algorithm to compute univocal extensions of omnitig hearts.
    type UnivocalExtensionStrategy: UnivocalExtensionAlgorithm<Graph, VecEdgeWalk<Graph>>;

    /// To a sequence of maximal non-trivial omnitigs add the maximal trivial omnitigs.
    /// The function should not compute any trivial omnitigs that are subwalks of maximal non-trivial omnitigs.
    fn compute_maximal_trivial_omnitigs(
        graph: &Graph,
        omnitigs: Omnitigs<Graph>,
    ) -> Omnitigs<Graph>;
}

/// The algorithm used to compute univocal extensions of omnitigs.
pub trait UnivocalExtensionAlgorithm<Graph: StaticGraph, ResultWalk: From<Vec<Graph::EdgeIndex>>> {
    /// Compute the univocal extension of a walk.
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::EdgeIndex]) -> ResultWalk;
}

/// A collection of node-centric omnitigs.
pub trait NodeCentricOmnitigs<
    Graph: GraphBase,
    NodeCentricOmnitigsSubsequence: for<'a> Sequence<'a, VecNodeWalk<Graph>, NodeCentricOmnitigsSubsequence> + ?Sized,
>:
    From<Vec<VecNodeWalk<Graph>>>
    + for<'a> Sequence<'a, VecNodeWalk<Graph>, NodeCentricOmnitigsSubsequence>
{
    /// Compute the trivial node-centric omnitigs in the given strongly connected graph.
    fn compute_trivial_node_centric_omnitigs(graph: &Graph) -> Self
    where
        Graph: StaticGraph,
    {
        DefaultTrivialNodeCentricOmnitigAlgorithm::<SccNodeCentricUnivocalExtensionStrategy>::compute_maximal_trivial_node_centric_omnitigs(graph, Vec::new()).into()
    }

    /// Compute the trivial node-centric omnitigs in the given graph that may not be strongly connected.
    fn compute_trivial_node_centric_omnitigs_non_scc(graph: &Graph) -> Self
    where
        Graph: StaticGraph,
    {
        DefaultTrivialNodeCentricOmnitigAlgorithm::<NonSccNodeCentricUnivocalExtensionStrategy>::compute_maximal_trivial_node_centric_omnitigs(graph, Vec::new()).into()
    }

    /// Retains only one direction of each pair of reverse-complemental omnitigs.
    ///
    /// Note: I am not sure if this method is correct in all cases, but it will panic if it finds a case where it is not correct.
    ///       For practical genomes it seems to work.
    fn remove_reverse_complements(&mut self, graph: &Graph)
    where
        Graph: StaticBigraph,
        Self: FromIterator<VecNodeWalk<Graph>>,
    {
        // Maps from nodes to omnitigs that start with this node.
        let mut first_node_map = vec![Vec::new(); graph.node_count()];
        for (i, omnitig) in self.iter().enumerate() {
            let first_node = omnitig.iter().next().expect("Omnitig is empty");
            first_node_map[first_node.as_usize()].push(i);
        }

        let mut retain_indices = Vec::with_capacity(self.len() / 2);
        for (i, omnitig) in self.iter().enumerate() {
            let reverse_complement_first_node = graph
                .mirror_node(*omnitig.last().expect("Omnitig is empty"))
                .expect("Node has no mirror");

            let mut reverse_complement_count = 0;
            let mut self_complemental = false;
            for &reverse_complement_candidate_index in
                &first_node_map[reverse_complement_first_node.as_usize()]
            {
                let reverse_complement_candidate = &self[reverse_complement_candidate_index];
                let is_reverse_complemental = omnitig
                    .iter()
                    .zip(reverse_complement_candidate.iter().rev())
                    .all(|(&n1, &n2)| n1 == graph.mirror_node(n2).expect("Node has no mirror"));

                if is_reverse_complemental {
                    debug_assert_eq!(omnitig.len(), reverse_complement_candidate.len(), "Walks are reverse complemental, but do not have the same length. This means one of them is not maximal.");
                    debug_assert_eq!(
                        reverse_complement_count, 0,
                        "Walk has more than one reverse complement."
                    );
                    reverse_complement_count += 1;

                    match reverse_complement_candidate_index.cmp(&i) {
                        Ordering::Less => retain_indices.push(reverse_complement_candidate_index),
                        Ordering::Equal => self_complemental = true,
                        Ordering::Greater => (),
                    }
                }

                if self_complemental {
                    retain_indices.push(reverse_complement_candidate_index);
                }
            }
        }

        retain_indices.sort_unstable();
        debug_assert!(
            retain_indices.windows(2).all(|w| w[0] < w[1]),
            "retain_indices contains duplicate walk"
        );

        let mut retained_omnitigs = Vec::new();
        for (i, omnitig) in self.iter().enumerate() {
            if retained_omnitigs.len() == retain_indices.len() {
                break;
            }

            if i == retain_indices[retained_omnitigs.len()] {
                retained_omnitigs.push(omnitig);
            }
        }

        *self = retained_omnitigs.into_iter().cloned().collect();
    }
}

impl<Graph: 'static + GraphBase> NodeCentricOmnitigs<Graph, [VecNodeWalk<Graph>]>
    for Vec<VecNodeWalk<Graph>>
{
}

/// A trait abstracting over the concrete algorithm used to compute maximal trivial node-centric omnitigs.
pub trait TrivialNodeCentricOmnitigAlgorithm<Graph: StaticGraph> {
    /// The algorithm to compute univocal extensions of node-centric omnitig hearts.
    type NodeCentricUnivocalExtensionStrategy: NodeCentricUnivocalExtensionAlgorithm<
        Graph,
        VecNodeWalk<Graph>,
    >;

    /// To a sequence of maximal non-trivial node-centric omnitigs add the maximal trivial node-centric omnitigs.
    /// The function should not compute any trivial node-centric omnitigs that are subwalks of maximal non-trivial node-centric omnitigs.
    fn compute_maximal_trivial_node_centric_omnitigs(
        graph: &Graph,
        omnitigs: Vec<VecNodeWalk<Graph>>,
    ) -> Vec<VecNodeWalk<Graph>>;
}

/// The algorithm used to compute univocal extensions of node-centric omnitigs.
pub trait NodeCentricUnivocalExtensionAlgorithm<
    Graph: StaticGraph,
    ResultWalk: From<Vec<Graph::NodeIndex>>,
>
{
    /// Compute the univocal extension of a node-centric walk.
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::NodeIndex]) -> ResultWalk;
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
        debug_assert_eq!(
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
