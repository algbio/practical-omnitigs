use crate::interface::DynamicGraph;
use rand::seq::IteratorRandom;
use rand::Rng;

/// Adds a binary tree to the given graph.
/// The first added node is the root of the tree.
/// A negative depth adds no nodes to the graph, a depth of 0 just the root, a depth of 1 the root an its children, and so on.
pub fn create_binary_tree<Graph: DynamicGraph>(
    graph: &mut Graph,
    depth: i32,
) -> Option<Graph::NodeIndex>
where
    Graph::NodeData: Default,
    Graph::EdgeData: Default,
{
    if depth < 0 {
        return None;
    }

    let root = graph.add_node(Default::default());
    create_binary_tree_recursively(graph, depth - 1, root);
    Some(root)
}

fn create_binary_tree_recursively<Graph: DynamicGraph>(
    graph: &mut Graph,
    depth: i32,
    root: Graph::NodeIndex,
) where
    Graph::NodeData: Default,
    Graph::EdgeData: Default,
{
    if depth < 0 {
        return;
    }

    let l = graph.add_node(Default::default());
    let r = graph.add_node(Default::default());
    graph.add_edge(root, l, Default::default());
    graph.add_edge(root, r, Default::default());
    create_binary_tree_recursively(graph, depth - 1, l);
    create_binary_tree_recursively(graph, depth - 1, r);
}

/// Computes the amount of edges in a graph with n nodes, given the hamiltonian edge factor c.
pub fn compute_m_from_n_and_c(n: usize, c: f64) -> usize {
    let node_amount_f64 = n as f64;
    let target_edge_amount =
        c * node_amount_f64 * (node_amount_f64.ln().max(1.0) + node_amount_f64.ln().ln().max(0.0));
    target_edge_amount.round() as usize
}

/// Creates a random hamiltonian graph with the given amount of nodes.
/// Assumes that the graph is empty.
/// The amount of arcs will be `c * n * (log(n) + log(log(n)))`, where `n` is the amount of nodes.
pub fn create_random_hamiltonian_graph<Graph: DynamicGraph, Random: Rng>(
    graph: &mut Graph,
    node_amount: usize,
    c: f64,
    random: &mut Random,
) where
    Graph::NodeData: Default,
    Graph::EdgeData: Default,
{
    if node_amount == 0 {
        return;
    }

    for _ in 0..node_amount {
        graph.add_node(Default::default());
    }
    for (n1, n2) in graph
        .node_indices()
        .take(graph.node_count() - 1)
        .zip(graph.node_indices().skip(1))
    {
        graph.add_edge(n1, n2, Default::default());
    }
    graph.add_edge(
        graph.node_indices().last().unwrap(),
        graph.node_indices().next().unwrap(),
        Default::default(),
    );

    let target_edge_amount = compute_m_from_n_and_c(node_amount, c);
    debug_assert!(
        target_edge_amount >= node_amount && target_edge_amount <= node_amount * (node_amount - 1),
        "node_amount <= target_edge_amount <= node_amount * (node_amount - 1): {} <= {} <= {} (c: {})",
        node_amount,
        target_edge_amount,
        node_amount * (node_amount - 1),
        c,
    );

    while graph.edge_count() < target_edge_amount {
        let n1 = graph.node_indices().choose(random).unwrap();
        let n2 = graph.node_indices().choose(random).unwrap();

        if n1 != n2 && !graph.contains_edge_between(n1, n2) {
            graph.add_edge(n1, n2, Default::default());
        }
    }
}

/// Creates a random graph with the given amount of nodes.
/// Assumes that the graph is empty.
/// The amount of arcs will be `c * n * (log(n) + log(log(n)))`, where `n` is the amount of nodes.
pub fn create_random_graph<Graph: DynamicGraph, Random: Rng>(
    graph: &mut Graph,
    node_amount: usize,
    c: f64,
    random: &mut Random,
) where
    Graph::NodeData: Default,
    Graph::EdgeData: Default,
{
    if node_amount == 0 {
        return;
    }

    for _ in 0..node_amount {
        graph.add_node(Default::default());
    }

    let target_edge_amount = compute_m_from_n_and_c(node_amount, c);
    debug_assert!(
        target_edge_amount >= node_amount && target_edge_amount <= node_amount * (node_amount - 1)
    );

    while graph.edge_count() < target_edge_amount {
        let n1 = graph.node_indices().choose(random).unwrap();
        let n2 = graph.node_indices().choose(random).unwrap();

        if n1 != n2 && !graph.contains_edge_between(n1, n2) {
            graph.add_edge(n1, n2, Default::default());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::create_binary_tree;
    use crate::implementation::petgraph_impl;
    use crate::interface::ImmutableGraphContainer;

    #[test]
    fn test_create_binary_tree_2() {
        let mut graph = petgraph_impl::new::<(), ()>();
        create_binary_tree(&mut graph, 2);
        debug_assert_eq!(graph.node_count(), 7);
        debug_assert_eq!(graph.edge_count(), 6);
    }
}
