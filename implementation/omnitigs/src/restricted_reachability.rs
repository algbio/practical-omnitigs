use traitgraph::implementation::incremental_subgraph::IncrementalSubgraph;
use traitgraph::interface::subgraph::DecoratingSubgraph;
use traitgraph::interface::{GraphBase, NodeOrEdge, StaticGraph};
use traitgraph_algo::traversal::{
    BackwardNeighborStrategy, BfsQueueStrategy, ForbiddenEdge, ForbiddenNode,
    ForwardNeighborStrategy, PreOrderTraversal, TraversalNeighborStrategy,
};

/// Returns the reachable subgraph from a node without using an edge.
pub fn compute_restricted_edge_reachability<
    'a,
    Graph: StaticGraph,
    NeighborStrategy: TraversalNeighborStrategy<'a, SubgraphType::ParentGraph>,
    SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
>(
    graph: SubgraphType::ParentGraphRef,
    start_node: <SubgraphType as GraphBase>::NodeIndex,
    forbidden_edge: <SubgraphType as GraphBase>::EdgeIndex,
) -> SubgraphType {
    let mut subgraph = SubgraphType::new_empty(graph);
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

/// Returns the reachable subgraph from a node without using an edge incrementally.
pub fn compute_incremental_restricted_forward_edge_reachability<'a, Graph: StaticGraph>(
    graph: &'a Graph,
    walk: &[Graph::EdgeIndex],
) -> IncrementalSubgraph<'a, Graph> {
    let mut subgraph = IncrementalSubgraph::new_with_incremental_steps(graph, walk.len());
    let mut traversal = PreOrderTraversal::<
        _,
        ForwardNeighborStrategy,
        BfsQueueStrategy,
        std::collections::VecDeque<_>,
    >::new_without_start(graph);

    let mut start_edge = *walk
        .first()
        .expect("Cannot compute hydrostructure from empty walk");
    for (edge_number, &edge) in walk.iter().enumerate().skip(1) {
        let start_node = graph.edge_endpoints(start_edge).to_node;
        traversal.continue_traversal_from(start_node);
        subgraph.set_current_step(edge_number);
        let forbidden_edge = ForbiddenEdge::new(edge);
        subgraph.add_edge(start_edge);
        start_edge = edge;

        while let Some(node_or_edge) = traversal.next_with_forbidden_subgraph(&forbidden_edge) {
            match node_or_edge {
                NodeOrEdge::Node(node) => subgraph.add_node(node),
                NodeOrEdge::Edge(edge) => subgraph.add_edge(edge),
            }
        }
    }

    subgraph
}

/// Returns the backwards reachable subgraph from a node without using an edge incrementally.
pub fn compute_incremental_restricted_backward_edge_reachability<'a, Graph: StaticGraph>(
    graph: &'a Graph,
    walk: &[Graph::EdgeIndex],
) -> IncrementalSubgraph<'a, Graph> {
    let mut subgraph = IncrementalSubgraph::new_with_incremental_steps(graph, walk.len());
    let mut traversal = PreOrderTraversal::<
        _,
        BackwardNeighborStrategy,
        BfsQueueStrategy,
        std::collections::VecDeque<_>,
    >::new_without_start(graph);

    let mut start_edge = *walk
        .last()
        .expect("Cannot compute hydrostructure from empty walk");
    for (edge_number, &edge) in walk.iter().rev().enumerate().skip(1) {
        let start_node = graph.edge_endpoints(start_edge).from_node;
        traversal.continue_traversal_from(start_node);
        subgraph.set_current_step(edge_number);
        let forbidden_edge = ForbiddenEdge::new(edge);
        subgraph.add_edge(start_edge);
        start_edge = edge;

        while let Some(node_or_edge) = traversal.next_with_forbidden_subgraph(&forbidden_edge) {
            match node_or_edge {
                NodeOrEdge::Node(node) => subgraph.add_node(node),
                NodeOrEdge::Edge(edge) => subgraph.add_edge(edge),
            }
        }
    }

    subgraph
}

/// Returns the reachable subgraph from a node without using a node.
pub fn compute_restricted_node_reachability<
    'a,
    Graph: StaticGraph,
    NeighborStrategy: TraversalNeighborStrategy<'a, SubgraphType::ParentGraph>,
    SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
>(
    graph: SubgraphType::ParentGraphRef,
    start_node: <SubgraphType as GraphBase>::NodeIndex,
    forbidden_node: <SubgraphType as GraphBase>::NodeIndex,
) -> SubgraphType {
    let mut subgraph = SubgraphType::new_empty(graph);
    let mut traversal = PreOrderTraversal::<
        _,
        NeighborStrategy,
        BfsQueueStrategy,
        std::collections::VecDeque<_>,
    >::new(graph, start_node);
    let forbidden_node = ForbiddenNode::new(forbidden_node);

    while let Some(node_or_edge) = traversal.next_with_forbidden_subgraph(&forbidden_node) {
        match node_or_edge {
            NodeOrEdge::Node(node) => subgraph.add_node(node),
            NodeOrEdge::Edge(edge) => subgraph.add_edge(edge),
        }
    }

    subgraph
}

/// Returns the forwards reachable subgraph from the tail of `edge` without using `edge`.
pub fn compute_restricted_forward_reachability<
    'a,
    Graph: StaticGraph,
    SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
>(
    graph: SubgraphType::ParentGraphRef,
    edge: <SubgraphType as GraphBase>::EdgeIndex,
) -> SubgraphType {
    let start_node = graph.edge_endpoints(edge).from_node;
    compute_restricted_edge_reachability::<_, ForwardNeighborStrategy, _>(graph, start_node, edge)
}

/// Returns the backwards reachable subgraph from the head of `edge` without using `edge`.
pub fn compute_restricted_backward_reachability<
    'a,
    Graph: StaticGraph,
    SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
>(
    graph: SubgraphType::ParentGraphRef,
    edge: <SubgraphType as GraphBase>::EdgeIndex,
) -> SubgraphType {
    let start_node = graph.edge_endpoints(edge).to_node;
    compute_restricted_edge_reachability::<_, BackwardNeighborStrategy, _>(graph, start_node, edge)
}

/// Returns the forwards reachable subgraph from `edge` without using the tail of `edge`.
pub fn compute_inverse_restricted_forward_reachability<
    'a,
    Graph: StaticGraph,
    SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
>(
    graph: SubgraphType::ParentGraphRef,
    edge: <SubgraphType as GraphBase>::EdgeIndex,
) -> SubgraphType {
    let forbidden_node = graph.edge_endpoints(edge).from_node;
    let start_node = graph.edge_endpoints(edge).to_node;

    // If the edge is a self loop.
    let mut result = if start_node == forbidden_node {
        SubgraphType::new_empty(graph)
    } else {
        compute_restricted_node_reachability::<_, ForwardNeighborStrategy, _>(
            graph,
            start_node,
            forbidden_node,
        )
    };

    result.add_edge(edge);
    result
}

/// Returns the backwards reachable subgraph from `edge` without using the head of `edge`.
pub fn compute_inverse_restricted_backward_reachability<
    'a,
    Graph: StaticGraph,
    SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
>(
    graph: SubgraphType::ParentGraphRef,
    edge: <SubgraphType as GraphBase>::EdgeIndex,
) -> SubgraphType {
    let forbidden_node = graph.edge_endpoints(edge).to_node;
    let start_node = graph.edge_endpoints(edge).from_node;

    // If the edge is a self loop.
    let mut result = if start_node == forbidden_node {
        SubgraphType::new_empty(graph)
    } else {
        compute_restricted_node_reachability::<_, BackwardNeighborStrategy, _>(
            graph,
            start_node,
            forbidden_node,
        )
    };

    result.add_edge(edge);
    result
}

/// Returns either the set of nodes and edges reachable from the first edge of aZb without using aZb as a subwalk,
/// or None, if the whole graph can be reached this way.
///
/// This computes `R⁺(aZb)` as defined in the hydrostructure paper.
/// If `Some` is returned, `aZb` is _bridge-like_, and otherwise it is _avertible_.
pub fn compute_hydrostructure_forward_reachability<
    'a,
    Graph: StaticGraph,
    SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
>(
    graph: SubgraphType::ParentGraphRef,
    azb: &[<<SubgraphType as DecoratingSubgraph>::ParentGraph as GraphBase>::EdgeIndex],
) -> Option<SubgraphType> {
    let a = *azb.iter().next().unwrap();
    let b = *azb.iter().last().unwrap();
    let start_node = graph.edge_endpoints(a).to_node;
    let mut subgraph =
        compute_restricted_edge_reachability::<_, ForwardNeighborStrategy, SubgraphType>(
            graph, start_node, b,
        );

    for &edge in azb.iter().take(azb.len() - 1) {
        let node = graph.edge_endpoints(edge).to_node;
        for incoming in graph.in_neighbors(node) {
            let incoming = incoming.edge_id;
            if incoming != edge && subgraph.contains_edge(incoming) {
                return None;
            }
        }
    }

    subgraph.add_edge(a);
    Some(subgraph)
}

/// Returns either the set of nodes and edges backwards reachable from the last edge of aZb without using aZb as a subwalk,
/// or None, if the whole graph can be reached this way.
///
/// This computes `R⁻(aZb)` as defined in the hydrostructure paper.
/// If `Some` is returned, `aZb` is _bridge-like_, and otherwise it is _avertible_.
pub fn compute_hydrostructure_backward_reachability<
    'a,
    Graph: StaticGraph,
    SubgraphType: DecoratingSubgraph<ParentGraph = Graph, ParentGraphRef = &'a Graph>,
>(
    graph: SubgraphType::ParentGraphRef,
    azb: &[<<SubgraphType as DecoratingSubgraph>::ParentGraph as GraphBase>::EdgeIndex],
) -> Option<SubgraphType> {
    let a = *azb.iter().next().unwrap();
    let b = *azb.iter().last().unwrap();
    let start_node = graph.edge_endpoints(b).from_node;
    let mut subgraph =
        compute_restricted_edge_reachability::<_, BackwardNeighborStrategy, SubgraphType>(
            graph, start_node, a,
        );

    for &edge in azb.iter().skip(1) {
        let node = graph.edge_endpoints(edge).from_node;
        for outgoing in graph.out_neighbors(node) {
            let outgoing = outgoing.edge_id;
            if outgoing != edge && subgraph.contains_edge(outgoing) {
                return None;
            }
        }
    }

    subgraph.add_edge(b);
    Some(subgraph)
}

#[cfg(test)]
mod tests {
    use crate::restricted_reachability::compute_restricted_backward_reachability;
    use crate::restricted_reachability::compute_restricted_forward_reachability;
    use traitgraph::implementation::bit_vector_subgraph::BitVectorSubgraph;
    use traitgraph::interface::subgraph::DecoratingSubgraph;
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
        let subgraph: BitVectorSubgraph<_> = compute_restricted_forward_reachability(&graph, e1);

        debug_assert_eq!(subgraph.node_count(), 2);
        debug_assert!(subgraph.contains_node(n0));
        debug_assert!(subgraph.contains_node(n2));

        debug_assert_eq!(subgraph.edge_count(), 2);
        debug_assert!(subgraph.contains_edge(e3));
        debug_assert!(subgraph.contains_edge(e5));
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
        let subgraph: BitVectorSubgraph<_> = compute_restricted_backward_reachability(&graph, e1);

        debug_assert_eq!(subgraph.node_count(), 2);
        debug_assert!(subgraph.contains_node(n0));
        debug_assert!(subgraph.contains_node(n2));

        debug_assert_eq!(subgraph.edge_count(), 2);
        debug_assert!(subgraph.contains_edge(e3));
        debug_assert!(subgraph.contains_edge(e5));
    }
}
