use bigraph::interface::static_bigraph::StaticEdgeCentricBigraph;
use bigraph::interface::BidirectedData;
use traitgraph::index::GraphIndex;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::{EdgeWalk, NodeWalk, VecEdgeWalk, VecNodeWalk};
use traitsequence::interface::Sequence;

/// Algorithms to count uncompacted unitigs in a graph.
pub mod uncompacted_unitigs;

/// A unitig stored as sequence of nodes.
#[derive(Clone)]
pub struct NodeUnitig<Graph: GraphBase> {
    /// Store the unitig as a node-centric walk.
    /// This allows unitigs that are a single node to be stored.
    walk: VecNodeWalk<Graph>,
    /// If a unitig has just a single edge, this edge might be a multiedge.
    /// In this case, this field is required to uniquely identify the unitig.
    single_edge_disambiguator: Option<Graph::EdgeIndex>,
}

impl<Graph: GraphBase> NodeUnitig<Graph> {
    /// Creates a new unitig from the given walk.
    /// Panics if the walk has length two nodes but no single edge is given, or if a single edge is given but the walk has a length other than two nodes.
    pub fn new(walk: VecNodeWalk<Graph>, single_edge: Option<Graph::EdgeIndex>) -> Self {
        debug_assert_eq!(walk.len() == 2, single_edge.is_some());
        Self {
            walk,
            single_edge_disambiguator: single_edge,
        }
    }

    /// Returns the single edge field of this unitig.
    /// The field is `Some` if this unitig has exactly two nodes.
    /// In this case, the two nodes might be connected by multiple edges, so then this field serves as disambiguation between the edges.
    pub fn single_edge(&self) -> &Option<Graph::EdgeIndex> {
        &self.single_edge_disambiguator
    }

    /// Extracts the node walk from this unitig, consuming the unitig.
    /// Note that the information for unitigs of length two with a multiedge is lost.
    pub fn into_node_walk(self) -> VecNodeWalk<Graph> {
        self.walk
    }
}

impl<'a, Graph: GraphBase> NodeWalk<'a, Graph, [Graph::NodeIndex]> for NodeUnitig<Graph> where
    Graph::NodeIndex: 'a
{
}

impl<'a, Graph: GraphBase> Sequence<'a, Graph::NodeIndex, [Graph::NodeIndex]> for NodeUnitig<Graph>
where
    Graph::NodeIndex: 'a,
{
    type Iterator = std::slice::Iter<'a, Graph::NodeIndex>;

    fn iter(&'a self) -> Self::Iterator {
        self.walk.iter()
    }

    fn len(&self) -> usize {
        self.walk.len()
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for NodeUnitig<Graph>
where
    VecNodeWalk<Graph>: std::ops::Index<IndexType>,
{
    type Output = <VecNodeWalk<Graph> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.walk.index(index)
    }
}

impl<Graph: GraphBase> std::fmt::Debug for NodeUnitig<Graph>
where
    Graph::NodeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "NodeUnitig[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

/// A unitig stored as sequence of edges.
#[derive(Clone)]
pub struct EdgeUnitig<Graph: GraphBase> {
    /// Store the unitig as an edge-centric walk.
    walk: VecEdgeWalk<Graph>,
}

impl<Graph: GraphBase> EdgeUnitig<Graph> {
    /// Creates a new unitig from the given walk.
    /// Panics if the walk is empty.
    pub fn new(walk: VecEdgeWalk<Graph>) -> Self {
        debug_assert!(!walk.is_empty());
        Self { walk }
    }
}

impl<'a, Graph: GraphBase> EdgeWalk<'a, Graph, [Graph::EdgeIndex]> for EdgeUnitig<Graph> where
    Graph::EdgeIndex: 'a
{
}

impl<'a, Graph: GraphBase> Sequence<'a, Graph::EdgeIndex, [Graph::EdgeIndex]> for EdgeUnitig<Graph>
where
    Graph::EdgeIndex: 'a,
{
    type Iterator = std::slice::Iter<'a, Graph::EdgeIndex>;

    fn iter(&'a self) -> Self::Iterator {
        self.walk.iter()
    }

    fn len(&self) -> usize {
        self.walk.len()
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for EdgeUnitig<Graph>
where
    VecEdgeWalk<Graph>: std::ops::Index<IndexType>,
{
    type Output = <VecEdgeWalk<Graph> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.walk.index(index)
    }
}

impl<Graph: GraphBase> std::fmt::Debug for EdgeUnitig<Graph>
where
    Graph::EdgeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "EdgeUnitig[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

/// A structure storing a vector of node-centric unitigs.
pub struct NodeUnitigs<Graph: GraphBase> {
    unitigs: Vec<NodeUnitig<Graph>>,
}

impl<Graph: StaticGraph> NodeUnitigs<Graph> {
    /// Computes the maximal unitigs of a graph.
    ///
    /// The unitigs are computed both node- and edge-centric.
    /// That means that parallel edges are counted as separate unitig each, and bivalent nodes are counted as unitigs as well.
    ///
    /// If the graph is a cycle, it outputs an arbitrary subwalk.
    pub fn compute(graph: &Graph) -> Self {
        let mut used_edges = vec![false; graph.edge_count()];
        let mut unitigs = Vec::new();

        // Add unitigs of length at least two nodes.
        for edge in graph.edge_indices() {
            if !used_edges[edge.as_usize()] {
                used_edges[edge.as_usize()] = true;

                let mut start_node = graph.edge_endpoints(edge).from_node;
                let mut end_node = graph.edge_endpoints(edge).to_node;
                let mut unitig = vec![end_node, start_node];

                while graph.is_biunivocal_node(start_node) && start_node != end_node {
                    let in_neighbor = graph.in_neighbors(start_node).next().unwrap();
                    start_node = in_neighbor.node_id;
                    used_edges[in_neighbor.edge_id.as_usize()] = true;
                    unitig.push(start_node);
                }

                unitig.reverse();

                while graph.is_biunivocal_node(end_node) && start_node != end_node {
                    let out_neighbor = graph.out_neighbors(end_node).next().unwrap();
                    end_node = out_neighbor.node_id;
                    used_edges[out_neighbor.edge_id.as_usize()] = true;
                    unitig.push(end_node);
                }

                let single_edge_disambiguator = if unitig.len() == 2 { Some(edge) } else { None };
                unitigs.push(NodeUnitig::new(unitig, single_edge_disambiguator));
            }
        }

        // Add single nodes
        for node in graph.node_indices() {
            if graph.is_bivalent_node(node) {
                unitigs.push(NodeUnitig::new(vec![node], None));
            }
        }

        Self { unitigs }
    }
}

/// A structure storing a vector of edge-centric unitigs.
pub struct EdgeUnitigs<Graph: GraphBase> {
    unitigs: Vec<EdgeUnitig<Graph>>,
}

impl<Graph: StaticGraph> EdgeUnitigs<Graph> {
    /// Computes the maximal unitigs of a graph.
    ///
    /// The unitigs are computed edge-centric.
    ///
    /// If the graph is a cycle, it outputs an arbitrary subwalk.
    pub fn compute(graph: &Graph) -> Self {
        let mut used_edges = vec![false; graph.edge_count()];
        let mut unitigs = Vec::new();

        // Add unitigs of length at least two nodes.
        for edge in graph.edge_indices() {
            if !used_edges[edge.as_usize()] {
                used_edges[edge.as_usize()] = true;

                let mut start_node = graph.edge_endpoints(edge).from_node;
                let mut end_node = graph.edge_endpoints(edge).to_node;
                let mut unitig = vec![edge];

                while graph.is_biunivocal_node(start_node) && start_node != end_node {
                    let in_neighbor = graph.in_neighbors(start_node).next().unwrap();
                    start_node = in_neighbor.node_id;
                    used_edges[in_neighbor.edge_id.as_usize()] = true;
                    unitig.push(in_neighbor.edge_id);
                }

                unitig.reverse();

                while graph.is_biunivocal_node(end_node) && start_node != end_node {
                    let out_neighbor = graph.out_neighbors(end_node).next().unwrap();
                    end_node = out_neighbor.node_id;
                    used_edges[out_neighbor.edge_id.as_usize()] = true;
                    unitig.push(out_neighbor.edge_id);
                }

                unitigs.push(EdgeUnitig::new(unitig));
            }
        }

        Self { unitigs }
    }

    /// Sorts the unitigs by length descending.
    pub fn sort_by_len_descending(&mut self) {
        self.unitigs.sort_by_key(|a| std::cmp::Reverse(a.len()));
    }
}

impl<Graph: StaticEdgeCentricBigraph> EdgeUnitigs<Graph>
where
    Graph::EdgeData: BidirectedData + Eq,
{
    /// Retains only one direction of each pair of reverse-complemental unitigs.
    pub fn remove_reverse_complements(&mut self, graph: &Graph) {
        // Maps from edges to unitigs that have this edge as first edge.
        let mut first_edge_map = vec![usize::max_value(); graph.edge_count()];
        for (i, unitig) in self.iter().enumerate() {
            let first_edge = unitig.iter().next().expect("Unitig is empty");
            debug_assert_eq!(
                first_edge_map[first_edge.as_usize()],
                usize::max_value(),
                "Found two unitigs starting with the same edge."
            );
            first_edge_map[first_edge.as_usize()] = i;
        }

        let mut retain_indices = Vec::with_capacity(self.len());
        for (i, unitig) in self.iter().enumerate() {
            let reverse_complement_first_edge = graph
                .mirror_edge_edge_centric(*unitig.iter().last().expect("Unitig his empty."))
                .expect("Edge has no reverse complement.");
            let reverse_complement_candidate_index =
                first_edge_map[reverse_complement_first_edge.as_usize()];
            if reverse_complement_candidate_index < i {
                let reverse_complement_candidate = &self[reverse_complement_candidate_index];
                for (edge, reverse_complement_edge) in
                    unitig.iter().zip(reverse_complement_candidate.iter().rev())
                {
                    debug_assert_eq!(
                        *edge,
                        graph
                            .mirror_edge_edge_centric(*reverse_complement_edge)
                            .expect("Edge has no reverse complement."),
                        "Found reverse complement candidate, but it is not a reverse complement."
                    );
                }
            } else {
                retain_indices.push(i);
            }
        }

        let mut unitigs = Vec::new();
        std::mem::swap(&mut unitigs, &mut self.unitigs);
        for (i, unitig) in unitigs.into_iter().enumerate() {
            if self.unitigs.len() == retain_indices.len() {
                break;
            }

            if i == retain_indices[self.unitigs.len()] {
                self.unitigs.push(unitig);
            }
        }
    }
}

impl<'a, Graph: 'a + GraphBase> Sequence<'a, NodeUnitig<Graph>, [NodeUnitig<Graph>]>
    for NodeUnitigs<Graph>
where
    Graph::NodeIndex: 'a,
{
    type Iterator = std::slice::Iter<'a, NodeUnitig<Graph>>;
    fn iter(&'a self) -> Self::Iterator {
        self.unitigs.iter()
    }

    fn len(&self) -> usize {
        self.unitigs.len()
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for NodeUnitigs<Graph>
where
    Vec<NodeUnitig<Graph>>: std::ops::Index<IndexType>,
{
    type Output = <Vec<NodeUnitig<Graph>> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.unitigs.index(index)
    }
}

impl<Graph: GraphBase> std::fmt::Debug for NodeUnitigs<Graph>
where
    Graph::NodeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "NodeUnitigs[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase> IntoIterator for NodeUnitigs<Graph> {
    type Item = NodeUnitig<Graph>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.unitigs.into_iter()
    }
}

impl<Graph: GraphBase> PartialEq for NodeUnitig<Graph>
where
    Graph::NodeIndex: PartialEq,
    Graph::EdgeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.walk == rhs.walk && self.single_edge_disambiguator == rhs.single_edge_disambiguator
    }
}

impl<Graph: GraphBase> Eq for NodeUnitig<Graph>
where
    Graph::NodeIndex: Eq,
    Graph::EdgeIndex: Eq,
{
}

impl<'a, Graph: 'a + GraphBase> Sequence<'a, EdgeUnitig<Graph>, [EdgeUnitig<Graph>]>
    for EdgeUnitigs<Graph>
where
    Graph::EdgeIndex: 'a,
{
    type Iterator = std::slice::Iter<'a, EdgeUnitig<Graph>>;

    fn iter(&'a self) -> Self::Iterator {
        self.unitigs.iter()
    }

    fn len(&self) -> usize {
        self.unitigs.len()
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for EdgeUnitigs<Graph>
where
    Vec<EdgeUnitig<Graph>>: std::ops::Index<IndexType>,
{
    type Output = <Vec<EdgeUnitig<Graph>> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.unitigs.index(index)
    }
}

impl<Graph: GraphBase> std::fmt::Debug for EdgeUnitigs<Graph>
where
    Graph::EdgeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "EdgeUnitigs[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase> IntoIterator for EdgeUnitigs<Graph> {
    type Item = EdgeUnitig<Graph>;
    type IntoIter = std::vec::IntoIter<EdgeUnitig<Graph>>;

    fn into_iter(self) -> Self::IntoIter {
        self.unitigs.into_iter()
    }
}

impl<Graph: GraphBase> PartialEq for EdgeUnitig<Graph>
where
    Graph::NodeIndex: PartialEq,
    Graph::EdgeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.walk == rhs.walk
    }
}

impl<Graph: GraphBase> Eq for EdgeUnitig<Graph>
where
    Graph::NodeIndex: Eq,
    Graph::EdgeIndex: Eq,
{
}

#[cfg(test)]
mod test {
    use super::{NodeUnitig, NodeUnitigs};
    use crate::unitigs::{EdgeUnitig, EdgeUnitigs};
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::{MutableGraphContainer, WalkableGraph};
    use traitsequence::interface::Sequence;

    #[test]
    fn test_node_unitig_computation() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let n6 = graph.add_node(6);
        let n7 = graph.add_node(7);
        let n8 = graph.add_node(8);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n3, n5, 14);
        graph.add_edge(n4, n8, 15);
        graph.add_edge(n5, n8, 16);
        let e7 = graph.add_edge(n8, n6, 17);
        let e8 = graph.add_edge(n8, n6, 175);
        graph.add_edge(n8, n7, 18);
        let e10 = graph.add_edge(n6, n0, 19);
        graph.add_edge(n7, n0, 20);

        let unitigs = NodeUnitigs::compute(&graph);
        let mut unitigs_iter = unitigs.iter();
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n0, n1, n2, n3]),
                None
            ))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n3, n4, n8]),
                None
            ))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n3, n5, n8]),
                None
            ))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n8, n6]),
                Some(e7)
            ))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n8, n6]),
                Some(e8)
            ))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n8, n7, n0]),
                None
            ))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n6, n0]),
                Some(e10)
            ))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(graph.create_node_walk(&[n8]), None))
        );
        debug_assert_eq!(unitigs_iter.next(), None);
    }

    #[test]
    fn test_edge_unitig_computation() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let n6 = graph.add_node(6);
        let n7 = graph.add_node(7);
        let n8 = graph.add_node(8);
        let e0 = graph.add_edge(n0, n1, 10);
        let e1 = graph.add_edge(n1, n2, 11);
        let e2 = graph.add_edge(n2, n3, 12);
        let e3 = graph.add_edge(n3, n4, 13);
        let e4 = graph.add_edge(n3, n5, 14);
        let e5 = graph.add_edge(n4, n8, 15);
        let e6 = graph.add_edge(n5, n8, 16);
        let e7 = graph.add_edge(n8, n6, 17);
        let e8 = graph.add_edge(n8, n6, 175);
        let e9 = graph.add_edge(n8, n7, 18);
        let e10 = graph.add_edge(n6, n0, 19);
        let e11 = graph.add_edge(n7, n0, 20);

        let unitigs = EdgeUnitigs::compute(&graph);
        let mut unitigs_iter = unitigs.iter();
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&EdgeUnitig::new(graph.create_edge_walk(&[e0, e1, e2])))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&EdgeUnitig::new(graph.create_edge_walk(&[e3, e5])))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&EdgeUnitig::new(graph.create_edge_walk(&[e4, e6])))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&EdgeUnitig::new(graph.create_edge_walk(&[e7])))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&EdgeUnitig::new(graph.create_edge_walk(&[e8])))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&EdgeUnitig::new(graph.create_edge_walk(&[e9, e11])))
        );
        debug_assert_eq!(
            unitigs_iter.next(),
            Some(&EdgeUnitig::new(graph.create_edge_walk(&[e10])))
        );
        debug_assert_eq!(unitigs_iter.next(), None);
    }
}
