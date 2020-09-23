use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecNodeWalk;

/// A macronode algorithm that requires the graph to be strongly connected.
pub mod strongly_connected_macronode_algorithm;

/// A struct containing the macronodes of an uncompressed graph, represented as walks through their uncompressed centers.
/// In an uncompressed graph, macronode centers are maximal unitigs with the property that their first node has outdegree = 1, and their last node has indegree = 1.
#[derive(Clone)]
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

    /// Returns the amount of macronodes in this struct.
    pub fn len(&self) -> usize {
        self.macronodes.len()
    }

    /// Returns true if there are no macronodes in this struct.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<'a, Graph: GraphBase> IntoIterator for &'a Macronodes<Graph> {
    type Item = &'a VecNodeWalk<Graph>;
    type IntoIter = std::slice::Iter<'a, VecNodeWalk<Graph>>;

    fn into_iter(self) -> Self::IntoIter {
        self.macronodes.iter()
    }
}

impl<Graph: GraphBase> std::fmt::Debug for Macronodes<Graph>
where
    Graph::NodeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Macronodes[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase> PartialEq for Macronodes<Graph>
where
    Graph::NodeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.macronodes == rhs.macronodes
    }
}

impl<Graph: GraphBase> Eq for Macronodes<Graph> where Graph::NodeIndex: Eq {}

/// A trait abstracting over the concrete algorithm used to compute macronodes.
pub trait MacronodeAlgorithm<Graph: StaticGraph> {
    /// Compute all macronodes in a graph.
    fn compute_macronodes(graph: &Graph) -> Macronodes<Graph>;
}
