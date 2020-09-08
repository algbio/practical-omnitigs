use crate::interface::DynamicGraph;

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

#[cfg(test)]
mod tests {
    use super::create_binary_tree;
    use crate::implementation::petgraph_impl;
    use crate::interface::ImmutableGraphContainer;

    #[test]
    fn test_create_binary_tree_2() {
        let mut graph = petgraph_impl::new::<(), ()>();
        create_binary_tree(&mut graph, 2);
        assert_eq!(graph.node_count(), 7);
        assert_eq!(graph.edge_count(), 6);
    }
}
