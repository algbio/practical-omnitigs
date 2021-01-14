use crate::algo::traversal::univocal_traversal::UnivocalIterator;
use crate::interface::{GraphBase, NodeOrEdge, StaticGraph};
use traitsequence::interface::Sequence;

/// A sequence of nodes in a graph, where each consecutive pair of nodes is connected by an edge.
pub trait NodeWalk<'a, Graph: GraphBase>: Sequence<'a, Graph::NodeIndex>
where
    Graph::NodeIndex: 'a,
{
    /// Computes the trivial heart of the walk, or returns `None` if the walk is non-trivial.
    /// The heart is returned as the index of the last split node and the index of the first join node of the walk.
    /// These are not part of the heart.
    fn compute_trivial_heart(&'a self, graph: &Graph) -> Option<(usize, usize)>
    where
        Graph: StaticGraph,
    {
        let mut last_split = 0;
        let mut first_join = self.len() - 1;
        for (i, &node) in self.iter().enumerate().take(self.len() - 1).skip(1) {
            if graph.is_split_node(node) {
                last_split = last_split.max(i);
            }
            if graph.is_join_node(node) {
                first_join = first_join.min(i);
            }
        }

        if last_split < first_join
            && (last_split > 0 || first_join < self.len() - 1 || self.len() == 2)
        {
            Some((last_split, first_join))
        } else {
            None
        }
    }

    /// Compute the amount of nodes in the trivial heart, or returns `None` if the walk is non-trivial.
    /// Recall that a heart is a walk from arc to arc.
    fn compute_trivial_heart_node_len(&'a self, graph: &Graph) -> Option<usize>
    where
        Graph: StaticGraph,
    {
        if let Some((start, end)) = self.compute_trivial_heart(graph) {
            Some(end - start - 1)
        } else {
            None
        }
    }

    /// Returns true if this walk is non-trivial.
    fn is_non_trivial(&'a self, graph: &Graph) -> bool
    where
        Graph: StaticGraph,
    {
        self.compute_trivial_heart(graph).is_none()
    }

    /// Computes the univocal extension of this walk.
    /// That is the concatenation LWR, where W is the walk, L the longest R-univocal walk to the first node of W and R the longest univocal walk from the last node of W.
    fn compute_univocal_extension<ResultWalk: From<Vec<Graph::NodeIndex>>>(
        &'a self,
        graph: &Graph,
    ) -> ResultWalk
    where
        Graph: StaticGraph,
    {
        self.compute_univocal_extension_with_original_offset(graph)
            .1
    }

    /// Computes the univocal extension of this walk.
    /// That is the concatenation LWR, where W is the walk, L the longest R-univocal walk to the first node of W and R the longest univocal walk from the last node of W.
    ///
    /// Additionally to the univocal extension, this function returns the offset of the original walk in the univocal extension as usize.
    fn compute_univocal_extension_with_original_offset<ResultWalk: From<Vec<Graph::NodeIndex>>>(
        &'a self,
        graph: &Graph,
    ) -> (usize, ResultWalk)
    where
        Graph: StaticGraph,
    {
        assert!(
            !self.is_empty(),
            "Cannot compute the univocal extension of an empty walk."
        );

        let mut result = Vec::new();
        for node_or_edge in UnivocalIterator::new_backward_without_start(
            graph,
            NodeOrEdge::Node(*self.first().unwrap()),
        ) {
            match node_or_edge {
                NodeOrEdge::Node(node) => result.push(node),
                NodeOrEdge::Edge(_) => {}
            }
        }

        result.reverse();
        let original_offset = result.len();
        result.extend(self.iter());

        for node_or_edge in UnivocalIterator::new_forward_without_start(
            graph,
            NodeOrEdge::Node(*self.last().unwrap()),
        ) {
            match node_or_edge {
                NodeOrEdge::Node(node) => result.push(node),
                NodeOrEdge::Edge(_) => {}
            }
        }

        (original_offset, ResultWalk::from(result))
    }

    /// Returns the edge walk represented by this node walk.
    /// If there is a consecutive pair of nodes with a multiedge, then None is returned.
    /// If this walk contains less than two nodes, then None is returned.
    /// If there is a consecutive pair of node not connected by an edge, then this method panics.
    fn clone_as_edge_walk<ResultWalk: From<Vec<Graph::EdgeIndex>>>(
        &'a self,
        graph: &Graph,
    ) -> Option<ResultWalk>
    where
        Graph: StaticGraph,
    {
        if self.len() < 2 {
            return None;
        }

        let mut walk = Vec::new();
        for node_pair in self.iter().take(self.len() - 1).zip(self.iter().skip(1)) {
            let from = *node_pair.0;
            let to = *node_pair.1;
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

        Some(ResultWalk::from(walk))
    }

    /// Returns true if this is a proper subwalk of the given walk.
    /// Proper means that the walks are not equal.
    fn is_proper_subwalk_of(&'a self, other: &Self) -> bool
    where
        Graph::NodeIndex: Eq,
    {
        self.is_proper_subsequence_of(other)
    }
}

/// A sequence of edges in a graph, where each consecutive pair of edges is connected by a node.
pub trait EdgeWalk<'a, Graph: GraphBase>: Sequence<'a, Graph::EdgeIndex>
where
    Graph::EdgeIndex: 'a,
{
    /// Computes the trivial heart of the walk, or returns `None` if the walk is non-trivial.
    /// The heart is returned as the index of the last split arc and the index of the first join arc of the walk.
    fn compute_trivial_heart(&'a self, graph: &Graph) -> Option<(usize, usize)>
    where
        Graph: StaticGraph,
    {
        let mut last_split = 0;
        let mut first_join = self.len() - 1;
        for (i, &edge) in self.iter().enumerate().take(self.len() - 1).skip(1) {
            if graph.is_split_edge(edge) {
                last_split = last_split.max(i);
            }
            if graph.is_join_edge(edge) {
                first_join = first_join.min(i);
            }
        }

        if last_split <= first_join
            && (last_split > 0 || first_join < self.len() - 1 || self.len() == 1)
        {
            Some((last_split, first_join))
        } else {
            None
        }
    }

    /// Compute the amount of edges in the trivial heart, or returns `None` if the walk is non-trivial.
    /// Recall that a heart is a walk from arc to arc.
    fn compute_trivial_heart_edge_len(&'a self, graph: &Graph) -> Option<usize>
    where
        Graph: StaticGraph,
    {
        if let Some((start, end)) = self.compute_trivial_heart(graph) {
            Some(end - start + 1)
        } else {
            None
        }
    }

    /// Returns true if this walk is non-trivial.
    fn is_non_trivial(&'a self, graph: &Graph) -> bool
    where
        Graph: StaticGraph,
    {
        self.compute_trivial_heart(graph).is_none()
    }

    /// Compute the univocal extension of a walk.
    /// That is the concatenation LWR, where W is the walk, L the longest R-univocal walk to the first edge of W and R the longest univocal walk from the last edge of W.
    fn compute_univocal_extension<ResultWalk: From<Vec<Graph::EdgeIndex>>>(
        &'a self,
        graph: &Graph,
    ) -> ResultWalk
    where
        Graph: StaticGraph,
    {
        self.compute_univocal_extension_with_original_offset(graph)
            .1
    }

    /// Compute the univocal extension of a walk.
    /// That is the concatenation LWR, where W is the walk, L the longest R-univocal walk to the first edge of W and R the longest univocal walk from the last edge of W.
    /// This variant handles not strongly connected graphs by forbidding L and R to repeat edges.
    fn compute_univocal_extension_non_scc<ResultWalk: From<Vec<Graph::EdgeIndex>>>(
        &'a self,
        graph: &Graph,
    ) -> ResultWalk
    where
        Graph: StaticGraph,
    {
        self.compute_univocal_extension_with_original_offset_non_scc(graph)
            .1
    }

    /// Compute the univocal extension of a walk.
    /// That is the concatenation LWR, where W is the walk, L the longest R-univocal walk to the first edge of W and R the longest univocal walk from the last edge of W.
    ///
    /// Additionally to the univocal extension, this function returns the offset of the original walk in the univocal extension as usize.
    fn compute_univocal_extension_with_original_offset<ResultWalk: From<Vec<Graph::EdgeIndex>>>(
        &'a self,
        graph: &Graph,
    ) -> (usize, ResultWalk)
    where
        Graph: StaticGraph,
    {
        assert!(
            !self.is_empty(),
            "Cannot compute the univocal extension of an empty walk."
        );

        let mut result = Vec::new();
        for node_or_edge in UnivocalIterator::new_backward_without_start(
            graph,
            NodeOrEdge::Edge(*self.first().unwrap()),
        ) {
            match node_or_edge {
                NodeOrEdge::Node(_) => {}
                NodeOrEdge::Edge(edge) => {
                    if &edge == self.first().unwrap() {
                        break;
                    } else {
                        result.push(edge)
                    }
                }
            }
        }

        result.reverse();
        let original_offset = result.len();
        result.extend(self.iter());

        for node_or_edge in UnivocalIterator::new_forward_without_start(
            graph,
            NodeOrEdge::Edge(*self.last().unwrap()),
        ) {
            match node_or_edge {
                NodeOrEdge::Node(_) => {}
                NodeOrEdge::Edge(edge) => {
                    if &edge == self.last().unwrap() {
                        break;
                    } else {
                        result.push(edge)
                    }
                }
            }
        }

        (original_offset, ResultWalk::from(result))
    }

    /// Compute the univocal extension of a walk.
    /// That is the concatenation LWR, where W is the walk, L the longest R-univocal walk to the first edge of W and R the longest univocal walk from the last edge of W.
    /// This variant handles not strongly connected graphs by forbidding L and R to repeat edges.
    ///
    /// Additionally to the univocal extension, this function returns the offset of the original walk in the univocal extension as usize.
    fn compute_univocal_extension_with_original_offset_non_scc<
        ResultWalk: From<Vec<Graph::EdgeIndex>>,
    >(
        &'a self,
        graph: &Graph,
    ) -> (usize, ResultWalk)
    where
        Graph: StaticGraph,
    {
        assert!(
            !self.is_empty(),
            "Cannot compute the univocal extension of an empty walk."
        );

        let mut result = Vec::new();
        for node_or_edge in UnivocalIterator::new_backward_without_start(
            graph,
            NodeOrEdge::Edge(*self.first().unwrap()),
        ) {
            match node_or_edge {
                NodeOrEdge::Node(_) => {}
                NodeOrEdge::Edge(edge) => {
                    if edge == *self.last().unwrap() || result.contains(&edge) {
                        break;
                    } else {
                        result.push(edge);
                    }
                }
            }
        }

        result.reverse();
        let original_offset = result.len();
        result.extend(self.iter());

        for node_or_edge in UnivocalIterator::new_forward_without_start(
            graph,
            NodeOrEdge::Edge(*self.last().unwrap()),
        ) {
            match node_or_edge {
                NodeOrEdge::Node(_) => {}
                NodeOrEdge::Edge(edge) => {
                    if edge == *self.first().unwrap()
                        || result[original_offset + self.len()..].contains(&edge)
                    {
                        break;
                    } else {
                        result.push(edge);
                    }
                }
            }
        }

        (original_offset, ResultWalk::from(result))
    }

    /// Returns the node walk represented by this edge walk.
    /// If this walk contains no edge, then None is returned.
    /// If there is a consecutive pair of edges not connected by a node, then this method panics.
    fn clone_as_node_walk<ResultWalk: From<Vec<Graph::NodeIndex>>>(
        &'a self,
        graph: &Graph,
    ) -> Option<ResultWalk>
    where
        Graph: StaticGraph,
    {
        if self.is_empty() {
            return None;
        }

        let mut walk = Vec::new();
        walk.push(
            graph
                .edge_endpoints(self.first().cloned().unwrap())
                .from_node,
        );
        for edge_pair in self.iter().take(self.len() - 1).zip(self.iter().skip(1)) {
            let node = graph.edge_endpoints(*edge_pair.0).to_node;
            assert_eq!(
                node,
                graph.edge_endpoints(*edge_pair.1).from_node,
                "Not a valid edge walk"
            );
            walk.push(node);
        }
        walk.push(graph.edge_endpoints(self.last().cloned().unwrap()).to_node);

        Some(ResultWalk::from(walk))
    }

    /// Returns true if this is a proper subwalk of the given walk.
    /// Proper means that the walks are not equal.
    fn is_proper_subwalk_of(&'a self, other: &Self) -> bool
    where
        Graph::EdgeIndex: Eq,
    {
        self.is_proper_subsequence_of(other)
    }
}

////////////////////
////// Slices //////
////////////////////

impl<'a, Graph: GraphBase> NodeWalk<'a, Graph> for [Graph::NodeIndex] where Graph::NodeIndex: 'a {}

impl<'a, Graph: GraphBase> EdgeWalk<'a, Graph> for [Graph::EdgeIndex] where Graph::EdgeIndex: 'a {}

/////////////////////////
////// VecNodeWalk //////
/////////////////////////

/// A node walk that is represented as a vector of node indices.
pub struct VecNodeWalk<Graph: GraphBase> {
    walk: Vec<Graph::NodeIndex>,
}

impl<Graph: GraphBase> VecNodeWalk<Graph> {
    /// Creates a new walk over the given node indices.
    pub fn new(walk: Vec<Graph::NodeIndex>) -> Self {
        Self { walk }
    }
}

impl<'a, Graph: GraphBase> NodeWalk<'a, Graph> for VecNodeWalk<Graph> where Graph::NodeIndex: 'a {}

impl<'a, Graph: GraphBase> Sequence<'a, Graph::NodeIndex> for VecNodeWalk<Graph>
where
    Graph::NodeIndex: 'a,
{
    type Iterator = std::slice::Iter<'a, Graph::NodeIndex>;
    type IteratorMut = std::slice::IterMut<'a, Graph::NodeIndex>;

    fn iter(&'a self) -> Self::Iterator {
        self.walk.iter()
    }

    fn iter_mut(&'a mut self) -> Self::IteratorMut {
        self.walk.iter_mut()
    }

    fn len(&self) -> usize {
        self.walk.len()
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

impl<Graph: GraphBase> Clone for VecNodeWalk<Graph>
where
    Vec<Graph::NodeIndex>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            walk: self.walk.clone(),
        }
    }
}

/////////////////////////
////// VecEdgeWalk //////
/////////////////////////

/// An edge walk that is represented as a vector of edge indices.
pub struct VecEdgeWalk<Graph: GraphBase> {
    walk: Vec<Graph::EdgeIndex>,
}

impl<Graph: GraphBase> VecEdgeWalk<Graph> {
    /// Creates a new walk over the given edge indices.
    pub fn new(walk: Vec<Graph::EdgeIndex>) -> Self {
        Self { walk }
    }
}

impl<'a, Graph: GraphBase> EdgeWalk<'a, Graph> for VecEdgeWalk<Graph> where Graph::EdgeIndex: 'a {}

impl<'a, Graph: GraphBase> Sequence<'a, Graph::EdgeIndex> for VecEdgeWalk<Graph>
where
    Graph::EdgeIndex: 'a,
{
    type Iterator = std::slice::Iter<'a, Graph::EdgeIndex>;
    type IteratorMut = std::slice::IterMut<'a, Graph::EdgeIndex>;

    fn iter(&'a self) -> Self::Iterator {
        self.walk.iter()
    }

    fn iter_mut(&'a mut self) -> Self::IteratorMut {
        self.walk.iter_mut()
    }

    fn len(&self) -> usize {
        self.walk.len()
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

impl<Graph: GraphBase> Clone for VecEdgeWalk<Graph>
where
    Vec<Graph::EdgeIndex>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            walk: self.walk.clone(),
        }
    }
}
