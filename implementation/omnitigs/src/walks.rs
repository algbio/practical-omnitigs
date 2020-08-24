use traitgraph::interface::GraphBase;

/// A sequence of nodes in a graph, where each consecutive pair of nodes is connected by an edge.
pub trait Walk<Graph: GraphBase> where for<'a> &'a Self: IntoIterator<Item = Graph::NodeIndex> {
}

/// A walk that is represented as a vector of node indices.
pub struct VecWalk<Graph: GraphBase> {
    walk: Vec<Graph::NodeIndex>,
}

impl<Graph: GraphBase> VecWalk<Graph> {
    pub fn new(walk: Vec<Graph::NodeIndex>) -> Self {
        Self {walk}
    }
}

impl<Graph: GraphBase> Walk<Graph> for VecWalk<Graph> {
}

impl<'a, Graph: GraphBase> IntoIterator for &'a VecWalk<Graph> {
    type Item = Graph::NodeIndex;
    type IntoIter = std::iter::Cloned<std::slice::Iter<'a, Self::Item>>;

    fn into_iter(self) -> Self::IntoIter {
        self.walk.iter().cloned()
    }
}

impl<Graph: GraphBase> From<Vec<Graph::NodeIndex>> for VecWalk<Graph> {
    fn from(vec: Vec<Graph::NodeIndex>) -> Self {
        Self::new(vec)
    }
}