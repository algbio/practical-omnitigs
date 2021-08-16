use crate::algo::traversal::{
    BackwardNeighborStrategy, ForwardNeighborStrategy, TraversalNeighborStrategy,
};
use crate::interface::{GraphBase, NavigableGraph, NodeOrEdge, StaticGraph};
use std::marker::PhantomData;

/// An iterator over the univocal extension of a node or edge.
/// The direction is defined by the `NeighborStrategy`.
///
/// Note that this iterator only works for directed graphs, as for undirected graphs it is unclear which node comes after an edge.
pub struct UnivocalIterator<'a, Graph: GraphBase, NeighborStrategy> {
    graph: &'a Graph,
    neighbor_strategy: PhantomData<NeighborStrategy>,
    current_node: Option<Graph::NodeIndex>,
    current_edge: Option<Graph::EdgeIndex>,
}

impl<
        'a,
        Graph: NavigableGraph<'a>,
        NeighborStrategy: 'a + TraversalNeighborStrategy<'a, Graph>,
    > UnivocalIterator<'a, Graph, NeighborStrategy>
{
    /// Create a new `UnivocalIterator` that iterates over the univocal extension of `start_element` including `start_element`, in the direction specified by `neighbor_strategy`.
    pub fn new(
        graph: &'a Graph,
        start_element: NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>,
    ) -> Self {
        let (current_node, current_edge) = match start_element {
            NodeOrEdge::Node(node) => (Some(node), None),
            NodeOrEdge::Edge(edge) => (None, Some(edge)),
        };
        Self {
            graph,
            neighbor_strategy: Default::default(),
            current_node,
            current_edge,
        }
    }

    /// Create a new `UnivocalIterator` that iterates over the univocal extension of `start_element` excluding `start_element`, in the direction specified by `neighbor_strategy`.
    pub fn new_without_start(
        graph: &'a Graph,
        start_element: NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>,
    ) -> impl 'a + Iterator<Item = NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>> {
        Self::new(graph, start_element).skip(1)
    }
}

impl<'a, Graph: NavigableGraph<'a>> UnivocalIterator<'a, Graph, ForwardNeighborStrategy>
where
    ForwardNeighborStrategy: TraversalNeighborStrategy<'a, Graph>,
{
    /// Create a new `UnivocalIterator` that iterates over the forward univocal extension of `start_element` including `start_element`.
    pub fn new_forward(
        graph: &'a Graph,
        start_element: NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>,
    ) -> Self {
        Self::new(graph, start_element)
    }

    /// Create a new `UnivocalIterator` that iterates over the forward univocal extension of `start_element` excluding `start_element`.
    pub fn new_forward_without_start(
        graph: &'a Graph,
        start_element: NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>,
    ) -> impl 'a + Iterator<Item = NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>> {
        Self::new_without_start(graph, start_element)
    }
}

impl<'a, Graph: NavigableGraph<'a>> UnivocalIterator<'a, Graph, BackwardNeighborStrategy>
where
    BackwardNeighborStrategy: TraversalNeighborStrategy<'a, Graph>,
{
    /// Create a new `UnivocalIterator` that iterates over the backward univocal extension of `start_element` including `start_element`.
    pub fn new_backward(
        graph: &'a Graph,
        start_element: NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>,
    ) -> Self {
        Self::new(graph, start_element)
    }

    /// Create a new `UnivocalIterator` that iterates over the backward univocal extension of `start_element` excluding `start_element`.
    pub fn new_backward_without_start(
        graph: &'a Graph,
        start_element: NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>,
    ) -> impl 'a + Iterator<Item = NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>> {
        Self::new_without_start(graph, start_element)
    }
}

impl<'a, Graph: NavigableGraph<'a>, NeighborStrategy: TraversalNeighborStrategy<'a, Graph>> Iterator
    for UnivocalIterator<'a, Graph, NeighborStrategy>
{
    type Item = NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.current_node, self.current_edge) {
            (Some(_), Some(edge)) => {
                self.current_edge = None;
                Some(NodeOrEdge::Edge(edge))
            }
            (Some(node), None) => {
                let mut next_neighbor_iter = NeighborStrategy::neighbor_iterator(self.graph, node);
                if let Some(neighbor) = next_neighbor_iter.next() {
                    let next_node = neighbor.node_id;
                    let next_edge = neighbor.edge_id;

                    if next_neighbor_iter.next().is_none() {
                        self.current_node = Some(next_node);
                        self.current_edge = Some(next_edge);
                    } else {
                        self.current_node = None;
                    }
                } else {
                    self.current_node = None;
                    self.current_edge = None;
                }

                Some(NodeOrEdge::Node(node))
            }
            (None, Some(edge)) => {
                let mut next_node_iter = NeighborStrategy::edge_neighbor_iterator(self.graph, edge);
                let next_node = next_node_iter
                    .next()
                    .expect("Edge does not have a node as successor.");
                debug_assert!(
                    next_node_iter.next().is_none(),
                    "Edge has more than one node as successor."
                );
                self.current_node = Some(next_node);
                self.current_edge = None;
                Some(NodeOrEdge::Edge(edge))
            }
            (None, None) => None,
        }
    }
}

/// Returns true if the given edge is self-bivalent in the given graph, i.e. its univocal extension repeats a node.
pub fn is_edge_self_bivalent<Graph: StaticGraph>(graph: &Graph, edge_id: Graph::EdgeIndex) -> bool {
    let forward_iter =
        UnivocalIterator::new_forward_without_start(graph, NodeOrEdge::Edge(edge_id));
    let mut backward_iter = UnivocalIterator::new_backward(graph, NodeOrEdge::Edge(edge_id));

    let mut last_element = None;
    for element in forward_iter {
        match element {
            NodeOrEdge::Edge(edge) => {
                if edge == edge_id {
                    return true;
                }
                last_element = Some(NodeOrEdge::Edge(edge));
            }
            node => last_element = Some(node),
        }
    }

    let last_element = last_element.expect(
        "Forward univocal extension is empty, but should at least contain the start edge itself",
    );
    backward_iter.any(|e| e == last_element)
}

#[cfg(test)]
mod tests {
    use crate::algo::traversal::univocal_traversal::is_edge_self_bivalent;
    use crate::implementation::petgraph_impl;
    use crate::interface::MutableGraphContainer;

    #[test]
    fn test_is_edge_self_bivalent_simple() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let n4 = graph.add_node(());
        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n0, n2, ());
        let e2 = graph.add_edge(n3, n0, ());
        let e3 = graph.add_edge(n4, n0, ());
        let e4 = graph.add_edge(n1, n3, ());
        let e5 = graph.add_edge(n1, n3, ());
        let e6 = graph.add_edge(n2, n4, ());
        let e7 = graph.add_edge(n2, n4, ());
        let e8 = graph.add_edge(n0, n0, ());

        debug_assert!(!is_edge_self_bivalent(&graph, e0));
        debug_assert!(!is_edge_self_bivalent(&graph, e1));
        debug_assert!(!is_edge_self_bivalent(&graph, e2));
        debug_assert!(!is_edge_self_bivalent(&graph, e3));
        debug_assert!(is_edge_self_bivalent(&graph, e4));
        debug_assert!(is_edge_self_bivalent(&graph, e5));
        debug_assert!(is_edge_self_bivalent(&graph, e6));
        debug_assert!(is_edge_self_bivalent(&graph, e7));
        debug_assert!(is_edge_self_bivalent(&graph, e8));
    }

    #[test]
    fn test_is_edge_self_bivalent_cycle() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n2, n0, ());

        debug_assert!(is_edge_self_bivalent(&graph, e0));
        debug_assert!(is_edge_self_bivalent(&graph, e1));
        debug_assert!(is_edge_self_bivalent(&graph, e2));
    }
}
