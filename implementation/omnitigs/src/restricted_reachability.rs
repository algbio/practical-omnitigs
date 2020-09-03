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
