use crate::interface::GraphBase;

/// A sequence of nodes in a graph, where each consecutive pair of nodes is connected by an edge.
pub trait NodeWalk<'a, Graph: GraphBase>: for<'b> From<&'b [Graph::NodeIndex]> {
    type Iter: Iterator<Item = Graph::NodeIndex>;

    fn iter(&'a self) -> Self::Iter;
}

/// A sequence of edges in a graph, where each consecutive pair of edges is connected by a node.
pub trait EdgeWalk<'a, Graph: GraphBase>: for<'b> From<&'b [Graph::EdgeIndex]> {
    type Iter: Iterator<Item = Graph::EdgeIndex>;

    fn iter(&'a self) -> Self::Iter;
}

/// A node walk that is represented as a vector of node indices.
#[derive(Clone, Debug)]
pub struct VecNodeWalk<Graph: GraphBase> {
    walk: Vec<Graph::NodeIndex>,
}

impl<Graph: GraphBase> VecNodeWalk<Graph> {
    pub fn new(walk: Vec<Graph::NodeIndex>) -> Self {
        Self { walk }
    }
}

impl<'a, Graph: GraphBase> NodeWalk<'a, Graph> for VecNodeWalk<Graph>
where
    Graph::NodeIndex: 'a,
{
    type Iter = std::iter::Cloned<std::slice::Iter<'a, Graph::NodeIndex>>;

    fn iter(&'a self) -> Self::Iter {
        self.walk.iter().cloned()
    }
}

impl<Graph: GraphBase> From<Vec<Graph::NodeIndex>> for VecNodeWalk<Graph> {
    fn from(vec: Vec<Graph::NodeIndex>) -> Self {
        Self::new(vec)
    }
}

impl<'a, Graph: GraphBase> From<&'a [Graph::NodeIndex]> for VecNodeWalk<Graph> {
    fn from(slice: &'a [Graph::NodeIndex]) -> Self {
        Self::new(slice.to_vec())
    }
}

impl<Graph: GraphBase> PartialEq for VecNodeWalk<Graph>
where
    Graph::NodeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.walk == rhs.walk
    }
}

impl<Graph: GraphBase> Eq for VecNodeWalk<Graph> where Graph::NodeIndex: Eq {}

/// An edge walk that is represented as a vector of edge indices.
#[derive(Clone, Debug)]
pub struct VecEdgeWalk<Graph: GraphBase> {
    walk: Vec<Graph::EdgeIndex>,
}

impl<Graph: GraphBase> VecEdgeWalk<Graph> {
    pub fn new(walk: Vec<Graph::EdgeIndex>) -> Self {
        Self { walk }
    }
}

impl<'a, Graph: GraphBase> EdgeWalk<'a, Graph> for VecEdgeWalk<Graph>
where
    Graph::EdgeIndex: 'a,
{
    type Iter = std::iter::Cloned<std::slice::Iter<'a, Graph::EdgeIndex>>;

    fn iter(&'a self) -> Self::Iter {
        self.walk.iter().cloned()
    }
}

impl<Graph: GraphBase> From<Vec<Graph::EdgeIndex>> for VecEdgeWalk<Graph> {
    fn from(vec: Vec<Graph::EdgeIndex>) -> Self {
        Self::new(vec)
    }
}

impl<'a, Graph: GraphBase> From<&'a [Graph::EdgeIndex]> for VecEdgeWalk<Graph> {
    fn from(slice: &'a [Graph::EdgeIndex]) -> Self {
        Self::new(slice.to_vec())
    }
}

impl<Graph: GraphBase> PartialEq for VecEdgeWalk<Graph>
where
    Graph::EdgeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.walk == rhs.walk
    }
}

impl<Graph: GraphBase> Eq for VecEdgeWalk<Graph> where Graph::EdgeIndex: Eq {}
