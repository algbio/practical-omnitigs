use traitgraph::algo::traversal::{
    BackwardNeighborStrategy, BfsQueueStrategy, ForbiddenEdge, ForwardNeighborStrategy,
    PreOrderTraversal, TraversalNeighborStrategy,
};
use traitgraph::implementation::bit_vector_subgraph::BitVectorSubgraph;
use traitgraph::interface::subgraph::Subgraph;
use traitgraph::interface::{NodeOrEdge, StaticGraph};

/// Returns the reachable subgraph from a node without using an edge e.
pub fn compute_restricted_reachability<
    'a,
    Graph: StaticGraph,
    NeighborStrategy: TraversalNeighborStrategy<'a, Graph>,
>(
    graph: &'a Graph,
    start_node: Graph::NodeIndex,
    forbidden_edge: Graph::EdgeIndex,
) -> impl Subgraph<Graph> {
    let mut subgraph = BitVectorSubgraph::new_empty(graph);
    let mut traversal = PreOrderTraversal::<
        _,
        NeighborStrategy,
        BfsQueueStrategy,
        std::collections::VecDeque<_>,
    >::new(graph, start_node);
    let forbidden_edge = ForbiddenEdge::new(forbidden_edge);

    while let Some(node_or_edge) = traversal.next_with_forbidden_subgraph(&forbidden_edge) {
        match node_or_edge {
            NodeOrEdge::Node(node) => subgraph.add_node(node),
            NodeOrEdge::Edge(edge) => subgraph.add_edge(edge),
        }
    }

    subgraph
}

/// Returns the forwards reachable subgraph from the tail of `edge` without using `edge`.
pub fn compute_restricted_forward_reachability<Graph: StaticGraph>(
    graph: &Graph,
    edge: Graph::EdgeIndex,
) -> impl Subgraph<Graph> {
    let start_node = graph.edge_endpoints(edge).from_node;
    compute_restricted_reachability::<_, ForwardNeighborStrategy>(graph, start_node, edge)
}

/// Returns the backwards reachable subgraph from the head of `edge` without using `edge`.
pub fn compute_restricted_backward_reachability<Graph: StaticGraph>(
    graph: &Graph,
    edge: Graph::EdgeIndex,
) -> impl Subgraph<Graph> {
    let start_node = graph.edge_endpoints(edge).to_node;
    compute_restricted_reachability::<_, BackwardNeighborStrategy>(graph, start_node, edge)
}

#[cfg(test)]
mod tests {
    use crate::restricted_reachability::compute_restricted_backward_reachability;
    use crate::restricted_reachability::compute_restricted_forward_reachability;
    use traitgraph::interface::subgraph::Subgraph;
    use traitgraph::interface::MutableGraphContainer;

    #[test]
    fn test_restricted_forward_reachability_simple() {
        let mut graph = traitgraph::implementation::petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let e1 = graph.add_edge(n0, n1, -1);
        let _e2 = graph.add_edge(n1, n1, -2);
        let e3 = graph.add_edge(n0, n0, -3);
        let _e4 = graph.add_edge(n1, n0, -4);
        let e5 = graph.add_edge(n0, n2, -5);
        let _e6 = graph.add_edge(n1, n2, -6);
        let subgraph = compute_restricted_forward_reachability(&graph, e1);

        assert_eq!(subgraph.node_count(), 2);
        assert!(subgraph.contains_node(n0));
        assert!(subgraph.contains_node(n2));

        assert_eq!(subgraph.edge_count(), 2);
        assert!(subgraph.contains_edge(e3));
        assert!(subgraph.contains_edge(e5));
    }

    #[test]
    fn test_restricted_backward_reachability_simple() {
        let mut graph = traitgraph::implementation::petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let e1 = graph.add_edge(n1, n0, -1);
        let _e2 = graph.add_edge(n1, n1, -2);
        let e3 = graph.add_edge(n0, n0, -3);
        let _e4 = graph.add_edge(n0, n1, -4);
        let e5 = graph.add_edge(n2, n0, -5);
        let _e6 = graph.add_edge(n2, n1, -6);
        let subgraph = compute_restricted_backward_reachability(&graph, e1);

        assert_eq!(subgraph.node_count(), 2);
        assert!(subgraph.contains_node(n0));
        assert!(subgraph.contains_node(n2));

        assert_eq!(subgraph.edge_count(), 2);
        assert!(subgraph.contains_edge(e3));
        assert!(subgraph.contains_edge(e5));
    }
}
