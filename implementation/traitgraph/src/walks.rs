use crate::interface::{GraphBase, StaticGraph};

/// A sequence of nodes in a graph, where each consecutive pair of nodes is connected by an edge.
pub trait NodeWalk<'a, Graph: GraphBase>:
    for<'b> From<&'b [Graph::NodeIndex]> + std::ops::Index<usize, Output = Graph::NodeIndex>
{
    /// The iterator type used to iterate over the nodes of this walk.
    type Iter: Iterator<Item = Graph::NodeIndex>;

    /// Returns an iterator over the nodes of this walk.
    fn iter(&'a self) -> Self::Iter;

    /// Returns the length of this walk as its amount of nodes.
    fn len(&'a self) -> usize {
        self.iter().count()
    }

    /// Returns true if this walk contains no nodes.
    fn is_empty(&'a self) -> bool {
        self.len() == 0
    }

    /// Returns the first node of this walk or `None`, if this walk is empty.
    fn first(&'a self) -> Option<Graph::NodeIndex> {
        self.iter().next()
    }

    /// Returns the last node of this walk or `None` if this walk is empty.
    fn last(&'a self) -> Option<Graph::NodeIndex> {
        self.iter().last()
    }
}

/// A sequence of edges in a graph, where each consecutive pair of edges is connected by a node.
pub trait EdgeWalk<'a, Graph: GraphBase>:
    for<'b> From<&'b [Graph::EdgeIndex]> + std::ops::Index<usize, Output = Graph::EdgeIndex>
{
    /// The iterator type used to iterate over the edges of this walk.
    type Iter: Iterator<Item = Graph::EdgeIndex>;

    /// Returns an iterator over the edges of this walk.
    fn iter(&'a self) -> Self::Iter;

    /// Returns the length of this walk as its amount of edges.
    fn len(&'a self) -> usize {
        self.iter().count()
    }

    /// Returns true if this walk contains no edges.
    fn is_empty(&'a self) -> bool {
        self.len() == 0
    }

    /// Returns the first edge of this walk or `None`, if this walk is empty.
    fn first(&'a self) -> Option<Graph::EdgeIndex> {
        self.iter().next()
    }

    /// Returns the last edge of this walk or `None` if this walk is empty.
    fn last(&'a self) -> Option<Graph::EdgeIndex> {
        self.iter().last()
    }
}

/////////////////////////
////// VecNodeWalk //////
/////////////////////////

/// A node walk that is represented as a vector of node indices.
#[derive(Clone)]
pub struct VecNodeWalk<Graph: GraphBase> {
    walk: Vec<Graph::NodeIndex>,
}

impl<Graph: GraphBase> VecNodeWalk<Graph> {
    /// Creates a new walk over the given node indices.
    pub fn new(walk: Vec<Graph::NodeIndex>) -> Self {
        Self { walk }
    }
}

impl<Graph: StaticGraph> VecNodeWalk<Graph> {
    /// Returns the edge walk represented by this node walk.
    /// If there is a consecutive pair of nodes with a multiedge, then None is returned.
    /// If this walk contains less than two nodes, then None is returned.
    /// If there is a consecutive pair of node not connected by an edge, then this method panics.
    pub fn clone_as_edge_walk<ResultWalk: for<'a> EdgeWalk<'a, Graph>>(
        &self,
        graph: &Graph,
    ) -> Option<ResultWalk> {
        if self.walk.len() < 2 {
            return None;
        }

        let mut walk = Vec::new();
        for node_pair in self.walk.windows(2) {
            let from = node_pair[0];
            let to = node_pair[1];
            let mut edges_between = graph.edges_between(from, to);

            if let Some(edge) = edges_between.next() {
                walk.push(edge);
            } else {
                panic!("Not a valid node walk");
            }
            if edges_between.next().is_some() {
                return None;
            }
        }

        Some(ResultWalk::from(&walk))
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

impl<Graph: GraphBase> std::fmt::Debug for VecNodeWalk<Graph>
where
    Graph::NodeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "NodeWalk[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for VecNodeWalk<Graph>
where
    Vec<Graph::NodeIndex>: std::ops::Index<IndexType>,
{
    type Output = <Vec<Graph::NodeIndex> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.walk.index(index)
    }
}

/////////////////////////
////// VecEdgeWalk //////
/////////////////////////

/// An edge walk that is represented as a vector of edge indices.
#[derive(Clone)]
pub struct VecEdgeWalk<Graph: GraphBase> {
    walk: Vec<Graph::EdgeIndex>,
}

impl<Graph: GraphBase> VecEdgeWalk<Graph> {
    /// Creates a new walk over the given edge indices.
    pub fn new(walk: Vec<Graph::EdgeIndex>) -> Self {
        Self { walk }
    }
}

impl<Graph: StaticGraph> VecEdgeWalk<Graph> {
    /// Returns the node walk represented by this edge walk.
    /// If this walk contains no edge, then None is returned.
    /// If there is a consecutive pair of edges not connected by a node, then this method panics.
    pub fn clone_as_node_walk<ResultWalk: for<'a> NodeWalk<'a, Graph>>(
        &self,
        graph: &Graph,
    ) -> Option<ResultWalk> {
        if self.walk.is_empty() {
            return None;
        }

        let mut walk = Vec::new();
        walk.push(
            graph
                .edge_endpoints(self.walk.first().cloned().unwrap())
                .from_node,
        );
        for edge_pair in self.walk.windows(2) {
            let node = graph.edge_endpoints(edge_pair[0]).to_node;
            assert_eq!(
                node,
                graph.edge_endpoints(edge_pair[1]).from_node,
                "Not a valid edge walk"
            );
            walk.push(node);
        }
        walk.push(
            graph
                .edge_endpoints(self.walk.last().cloned().unwrap())
                .to_node,
        );

        Some(ResultWalk::from(&walk))
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

impl<Graph: GraphBase> std::fmt::Debug for VecEdgeWalk<Graph>
where
    Graph::EdgeIndex: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "EdgeWalk[")?;
        if let Some(first) = self.iter().next() {
            write!(f, "{:?}", first)?;
        }
        for edge in self.iter().skip(1) {
            write!(f, ", {:?}", edge)?;
        }
        write!(f, "]")
    }
}

impl<Graph: GraphBase, IndexType> std::ops::Index<IndexType> for VecEdgeWalk<Graph>
where
    Vec<Graph::EdgeIndex>: std::ops::Index<IndexType>,
{
    type Output = <Vec<Graph::EdgeIndex> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.walk.index(index)
    }
}
