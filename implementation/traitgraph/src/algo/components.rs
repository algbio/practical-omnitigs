use super::traversal::{
    PostOrderForwardDfs, PreOrderBackwardBfs, PreOrderForwardBfs, PreOrderUndirectedBfs,
};
use crate::algo::traversal::ForbiddenEdge;
use crate::index::GraphIndex;
use crate::index::OptionalGraphIndex;
use crate::interface::NodeOrEdge;
use crate::interface::{DynamicGraph, MutableGraphContainer, StaticGraph};
use std::collections::LinkedList;

/// Returns the weakly connected components of a graph.
///
/// If the graph is empty, no WCCs are returned.
/// Otherwise, the WCCs are cloned into new graphs, without preserving the node or edge indices.
pub fn decompose_weakly_connected_components<Graph: Default + DynamicGraph>(
    graph: &Graph,
) -> Vec<Graph>
where
    Graph::NodeData: Clone,
    Graph::EdgeData: Clone,
{
    let mut result = Vec::new();
    let mut nodes: Vec<_> = graph.node_indices().collect();
    let mut visited = vec![false; graph.node_count()];
    let mut bfs = PreOrderUndirectedBfs::new_without_start(graph);

    while !nodes.is_empty() {
        let start = nodes.pop().unwrap();
        if visited[start.as_usize()] {
            continue;
        }

        let rank_offset = bfs.continue_traversal_from(start).as_usize();
        let mut subgraph = Graph::default();

        while let Some(node) = bfs.next() {
            let node = if let NodeOrEdge::Node(node) = node {
                node
            } else {
                continue;
            };
            visited[node.as_usize()] = true;
            //println!("add_node: {:?}", node);
            let subnode = subgraph.add_node(graph.node_data(node).clone());

            for out_neighbor in graph.out_neighbors(node) {
                let neighbor_id = out_neighbor.node_id;
                debug_assert!(
                    graph.contains_edge_between(node, neighbor_id),
                    "f: Edge missing: ({:?}, {:?})",
                    node,
                    neighbor_id
                );
                let edge_id = out_neighbor.edge_id;
                if let Some(subneighbor) = bfs.rank_of(neighbor_id) {
                    let subneighbor = (subneighbor.as_usize() - rank_offset).into();
                    if subgraph.contains_node_index(subneighbor) {
                        let edge_data = graph.edge_data(edge_id).clone();
                        //println!("f: ({:?}, {:?}) becomes ({:?}, {:?})", node, neighbor_id, subnode, subneighbor);
                        subgraph.add_edge(subnode, subneighbor, edge_data);
                    }
                }
            }
            for in_neighbor in graph.in_neighbors(node) {
                let neighbor_id = in_neighbor.node_id;

                // Handle self loops only once in the forward case.
                // Otherwise, two identical versions of the loop would be added to the subgraph.
                if neighbor_id == node {
                    continue;
                }
                debug_assert!(
                    graph.contains_edge_between(neighbor_id, node),
                    "r: Edge missing: ({:?}, {:?})",
                    neighbor_id,
                    node
                );
                let edge_id = in_neighbor.edge_id;

                if let Some(subneighbor) = bfs.rank_of(neighbor_id) {
                    let subneighbor = (subneighbor.as_usize() - rank_offset).into();
                    if subgraph.contains_node_index(subneighbor) {
                        let edge_data = graph.edge_data(edge_id).clone();
                        //println!("r: ({:?}, {:?}) becomes ({:?}, {:?})", neighbor_id, node, subneighbor, subnode);
                        subgraph.add_edge(subneighbor, subnode, edge_data);
                    }
                }
            }
        }

        result.push(subgraph);
    }

    result
}

/// Returns true if the graph is strongly connected.
pub fn is_strongly_connected<Graph: StaticGraph>(graph: &Graph) -> bool {
    if graph.is_empty() {
        return true;
    }

    let traversal = PreOrderForwardBfs::new(graph, graph.node_indices().next().unwrap());
    let mut traversal_node_count = 0;

    for node_or_edge in traversal {
        if let NodeOrEdge::Node(_) = node_or_edge {
            traversal_node_count += 1;
        }
    }

    if traversal_node_count != graph.node_count() {
        debug_assert!(traversal_node_count < graph.node_count());
        return false;
    }

    let traversal = PreOrderBackwardBfs::new(graph, graph.node_indices().next().unwrap());
    let mut traversal_node_count = 0;

    for node_or_edge in traversal {
        if let NodeOrEdge::Node(_) = node_or_edge {
            traversal_node_count += 1;
        }
    }

    if traversal_node_count != graph.node_count() {
        debug_assert!(traversal_node_count < graph.node_count());
        return false;
    }

    true
}

/// Returns true if the given edge is a strong bridge.
/// Note that this function assumes that the graph is strongly connected, and always returns true otherwise.
pub fn is_strong_bridge<Graph: StaticGraph>(graph: &Graph, edge: Graph::EdgeIndex) -> bool {
    debug_assert!(is_strongly_connected(graph));

    let mut traversal = PreOrderForwardBfs::new(graph, graph.node_indices().next().unwrap());
    let forbidden_edge = ForbiddenEdge::new(edge);
    let mut traversal_node_count = 0;

    while let Some(node_or_edge) = traversal.next_with_forbidden_subgraph(&forbidden_edge) {
        if let NodeOrEdge::Node(_) = node_or_edge {
            traversal_node_count += 1;
        }
    }

    if traversal_node_count != graph.node_count() {
        debug_assert!(traversal_node_count < graph.node_count());
        return true;
    }

    let mut traversal = PreOrderBackwardBfs::new(graph, graph.node_indices().next().unwrap());
    let mut traversal_node_count = 0;

    while let Some(node_or_edge) = traversal.next_with_forbidden_subgraph(&forbidden_edge) {
        if let NodeOrEdge::Node(_) = node_or_edge {
            traversal_node_count += 1;
        }
    }

    if traversal_node_count != graph.node_count() {
        debug_assert!(traversal_node_count < graph.node_count());
        return true;
    }

    false
}

/// Returns the strongly connected components of a graph.
///
/// If the graph is empty, no SCCs are returned.
/// Otherwise, an array is returned that maps each node to a root node representing its SCC.
/// Node that if the node ids are not consecutive, this mapping is still returned as consecutive array.
pub fn decompose_strongly_connected_components<Graph: StaticGraph>(
    graph: &Graph,
) -> Vec<Graph::NodeIndex> {
    let mut result: Vec<_> = graph.node_indices().collect();
    let mut nodes = LinkedList::new();
    let mut visited = vec![false; graph.node_count()];
    // 0 will be overridden with the first reset.
    let mut dfs = PostOrderForwardDfs::new_without_start(graph);

    for node in graph.node_indices() {
        if !visited[node.as_usize()] {
            dfs.continue_traversal_from(node);

            while let Some(node) = dfs.next(graph) {
                visited[node.as_usize()] = true;
                nodes.push_front(node);
            }
        }
    }

    //println!("nodes: {:?}", nodes);

    let mut bfs = PreOrderBackwardBfs::new_without_start(graph);
    for root_node in nodes {
        if visited[root_node.as_usize()] {
            //println!("Reverse processing {:?}", root_node);
            bfs.continue_traversal_from(root_node);

            for node in &mut bfs {
                let node = if let NodeOrEdge::Node(node) = node {
                    node
                } else {
                    continue;
                };
                visited[node.as_usize()] = false;
                result[node.as_usize()] = root_node;
            }
        }
    }

    //println!("result: {:?}", result);
    result
}

/// Extract the subgraphs of the given graph according to the given node_mapping.
///
/// The node indices of the graph are assumed to match the indices of the vector given as node mapping.
/// The return value is a vector of graphs of which each is the induced subgraph of a set of nodes with the same mapped value.
pub fn extract_subgraphs_from_node_mapping<Graph: Default + MutableGraphContainer + StaticGraph>(
    graph: &Graph,
    node_mapping: &[Graph::NodeIndex],
) -> Vec<Graph>
where
    Graph::NodeData: Clone,
    Graph::EdgeData: Clone,
{
    let mut result = Vec::new();
    let mut extracted_nodes = vec![Graph::OptionalNodeIndex::new_none(); graph.node_count()];
    let mut id_map = Vec::new();

    for node in graph.node_indices() {
        if !extracted_nodes[node.as_usize()].is_valid() {
            //println!("Processing {:?}", node);

            let root_node = node_mapping[node.as_usize()];
            let mut subgraph = Graph::default();
            id_map.clear();

            for node in graph
                .node_indices()
                .skip(node.as_usize())
                .filter(|n| node_mapping[n.as_usize()] == root_node)
            {
                //println!("Adding node {:?}", node);
                let subgraph_node = subgraph.add_node(graph.node_data(node).clone());
                extracted_nodes[node.as_usize()] = subgraph_node.into();
                id_map.push(node);
            }

            for subgraph_node in subgraph.node_indices() {
                let node = id_map[subgraph_node.as_usize()];
                for neighbor in graph.out_neighbors(node) {
                    let edge = neighbor.edge_id;
                    let neighbor_node = neighbor.node_id;

                    if node_mapping[neighbor_node.as_usize()] == root_node {
                        let subgraph_neighbor_node =
                            extracted_nodes[neighbor_node.as_usize()].into().unwrap();
                        let edge_data = graph.edge_data(edge);
                        subgraph.add_edge(subgraph_node, subgraph_neighbor_node, edge_data.clone());
                    }
                }
            }

            result.push(subgraph);
        }
    }

    result
}

/// Returns true if the given graph is a cycle.
pub fn is_cycle<Graph: StaticGraph>(graph: &Graph) -> bool {
    is_strongly_connected(graph) && graph.node_count() == graph.edge_count()
}

#[cfg(test)]
mod tests {
    use crate::algo::components::{
        decompose_strongly_connected_components, decompose_weakly_connected_components,
        extract_subgraphs_from_node_mapping, is_strongly_connected,
    };
    use crate::implementation::petgraph_impl;
    use crate::interface::{ImmutableGraphContainer, MutableGraphContainer};

    #[test]
    fn test_decompose_weakly_connected_components_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let _e0 = graph.add_edge(n0, n1, 10);
        let _e1 = graph.add_edge(n1, n2, 11);
        let _e15 = graph.add_edge(n1, n2, 115);
        let _e2 = graph.add_edge(n2, n3, 12);
        let _e3 = graph.add_edge(n3, n4, 13);
        let _e4 = graph.add_edge(n4, n0, 14);
        let _e45 = graph.add_edge(n4, n0, 145);
        let _e5 = graph.add_edge(n1, n0, 15);
        let _e6 = graph.add_edge(n2, n1, 16);
        let _e7 = graph.add_edge(n3, n2, 17);
        let _e8 = graph.add_edge(n4, n3, 18);
        let _e9 = graph.add_edge(n0, n4, 19);
        let _e10 = graph.add_edge(n2, n2, 20);
        let result = decompose_weakly_connected_components(&graph);
        debug_assert_eq!(result.len(), 1);
        let result = result.first().unwrap();
        debug_assert_eq!(result.node_count(), graph.node_count());
        debug_assert_eq!(result.edge_count(), graph.edge_count());

        let mut node_data: Vec<_> = graph
            .node_indices()
            .map(|i| result.node_data(i))
            .copied()
            .collect();
        node_data.sort_unstable();
        debug_assert_eq!(node_data, vec![0, 1, 2, 3, 4]);

        let mut edge_data: Vec<_> = graph
            .edge_indices()
            .map(|i| result.edge_data(i))
            .copied()
            .collect();
        edge_data.sort_unstable();
        debug_assert_eq!(
            edge_data,
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 115, 145]
        );
    }

    #[test]
    fn test_decompose_weakly_connected_components_one_wc_with_two_sccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let _e0 = graph.add_edge(n0, n1, 10);
        let _e1 = graph.add_edge(n1, n2, 11);
        let _e2 = graph.add_edge(n2, n3, 12);
        let _e3 = graph.add_edge(n3, n4, 13);
        let _e4 = graph.add_edge(n1, n0, 15);
        let _e5 = graph.add_edge(n3, n2, 17);
        let _e6 = graph.add_edge(n4, n3, 18);
        let _e7 = graph.add_edge(n2, n2, 20);
        let _e8 = graph.add_edge(n2, n2, 21);
        let result = decompose_weakly_connected_components(&graph);
        debug_assert_eq!(result.len(), 1);
        let result = result.first().unwrap();
        debug_assert_eq!(result.node_count(), graph.node_count());
        debug_assert_eq!(result.edge_count(), graph.edge_count());

        let mut node_data: Vec<_> = graph
            .node_indices()
            .map(|i| result.node_data(i))
            .copied()
            .collect();
        node_data.sort_unstable();
        debug_assert_eq!(node_data, vec![0, 1, 2, 3, 4]);

        let mut edge_data: Vec<_> = graph
            .edge_indices()
            .map(|i| result.edge_data(i))
            .copied()
            .collect();
        edge_data.sort_unstable();
        debug_assert_eq!(edge_data, vec![10, 11, 12, 13, 15, 17, 18, 20, 21]);
    }

    #[test]
    fn test_decompose_weakly_connected_components_multiple_wccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        graph.add_edge(n0, n0, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n1, n2, 115);
        graph.add_edge(n3, n4, 12);
        graph.add_edge(n4, n5, 13);
        let result = decompose_weakly_connected_components(&graph);
        debug_assert_eq!(result.len(), 3);
        let first = result[2].clone();
        let second = result[1].clone();
        let third = result[0].clone();
        debug_assert_eq!(first.node_count(), 1);
        debug_assert_eq!(first.edge_count(), 1);
        debug_assert_eq!(second.node_count(), 2);
        debug_assert_eq!(second.edge_count(), 2);
        debug_assert_eq!(third.node_count(), 3);
        debug_assert_eq!(third.edge_count(), 2);

        debug_assert_eq!(first.node_data(0.into()), &0);
        debug_assert_eq!(second.node_data(1.into()), &1);
        debug_assert_eq!(second.node_data(0.into()), &2);
        debug_assert_eq!(third.node_data(2.into()), &3);
        debug_assert_eq!(third.node_data(1.into()), &4);
        debug_assert_eq!(third.node_data(0.into()), &5);

        debug_assert_eq!(first.edge_data(0.into()), &10);
        debug_assert_eq!(second.edge_data(0.into()), &115);
        debug_assert_eq!(second.edge_data(1.into()), &11);
        debug_assert_eq!(third.edge_data(1.into()), &12);
        debug_assert_eq!(third.edge_data(0.into()), &13);
    }

    #[test]
    fn test_decompose_weakly_connected_components_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        let result = decompose_weakly_connected_components(&graph);
        debug_assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_decompose_weakly_connected_components_nearly_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let _e0 = graph.add_edge(n0, n1, 10);
        let _e05 = graph.add_edge(n0, n1, 105);
        let _e1 = graph.add_edge(n1, n2, 11);
        let _e2 = graph.add_edge(n2, n3, 12);
        let _e3 = graph.add_edge(n3, n4, 13);
        let _e4 = graph.add_edge(n4, n1, 14);
        let _e5 = graph.add_edge(n1, n4, 15);
        let _e6 = graph.add_edge(n2, n1, 16);
        let _e7 = graph.add_edge(n3, n2, 17);
        let _e8 = graph.add_edge(n4, n3, 18);
        let _e9 = graph.add_edge(n0, n4, 19);
        let _e10 = graph.add_edge(n2, n2, 20);
        let result = decompose_weakly_connected_components(&graph);
        debug_assert_eq!(result.len(), 1);
        let result = result.first().unwrap();
        debug_assert_eq!(result.node_count(), graph.node_count());
        debug_assert_eq!(result.edge_count(), graph.edge_count());

        let mut node_data: Vec<_> = graph
            .node_indices()
            .map(|i| result.node_data(i))
            .copied()
            .collect();
        node_data.sort_unstable();
        debug_assert_eq!(node_data, vec![0, 1, 2, 3, 4]);

        let mut edge_data: Vec<_> = graph
            .edge_indices()
            .map(|i| result.edge_data(i))
            .copied()
            .collect();
        edge_data.sort_unstable();
        debug_assert_eq!(
            edge_data,
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 105]
        );
    }

    #[test]
    fn test_decompose_weakly_connected_components_nearly_scc_reverse() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let _e0 = graph.add_edge(n4, n1, 10);
        let _e1 = graph.add_edge(n1, n2, 11);
        let _e2 = graph.add_edge(n2, n3, 12);
        let _e3 = graph.add_edge(n3, n4, 13);
        let _e4 = graph.add_edge(n4, n0, 14);
        let _e5 = graph.add_edge(n1, n0, 15);
        let _e55 = graph.add_edge(n1, n0, 155);
        let _e6 = graph.add_edge(n2, n1, 16);
        let _e7 = graph.add_edge(n3, n2, 17);
        let _e8 = graph.add_edge(n4, n3, 18);
        let _e9 = graph.add_edge(n1, n4, 19);
        let _e10 = graph.add_edge(n2, n2, 20);
        let result = decompose_weakly_connected_components(&graph);
        debug_assert_eq!(result.len(), 1);
        let result = result.first().unwrap();
        debug_assert_eq!(result.node_count(), graph.node_count());
        debug_assert_eq!(result.edge_count(), graph.edge_count());

        let mut node_data: Vec<_> = graph
            .node_indices()
            .map(|i| result.node_data(i))
            .copied()
            .collect();
        node_data.sort_unstable();
        debug_assert_eq!(node_data, vec![0, 1, 2, 3, 4]);

        let mut edge_data: Vec<_> = graph
            .edge_indices()
            .map(|i| result.edge_data(i))
            .copied()
            .collect();
        edge_data.sort_unstable();
        debug_assert_eq!(
            edge_data,
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 155]
        );
    }

    #[test]
    fn test_scc_check_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n0, 14);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n0, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(is_strongly_connected(&graph));
    }

    #[test]
    fn test_scc_check_one_wc_with_two_sccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
    }

    #[test]
    fn test_scc_check_multiple_wccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        graph.add_edge(n0, n0, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n1, 13);
        graph.add_edge(n3, n4, 12);
        graph.add_edge(n3, n4, 125);
        graph.add_edge(n4, n5, 13);
        graph.add_edge(n5, n3, 13);
        debug_assert!(!is_strongly_connected(&graph));
    }

    #[test]
    fn test_scc_check_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        debug_assert!(is_strongly_connected(&graph));
    }

    #[test]
    fn test_scc_check_nearly_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n1, 14);
        graph.add_edge(n1, n4, 15);
        graph.add_edge(n1, n4, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n0, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
    }

    #[test]
    fn test_scc_check_nearly_scc_reverse() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n4, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n0, 14);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n1, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
    }

    /////////////////////////////////////////////////////////

    #[test]
    fn test_decompose_sccs_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n0, 14);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n0, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(is_strongly_connected(&graph));
        debug_assert_eq!(decompose_strongly_connected_components(&graph), vec![n0; 5]);
    }

    #[test]
    fn test_decompose_sccs_one_wc_with_two_sccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 1);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
        debug_assert_eq!(
            decompose_strongly_connected_components(&graph),
            vec![n0, n0, n2, n2, n2]
        );
    }

    #[test]
    fn test_decompose_sccs_multiple_wccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        graph.add_edge(n0, n0, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n1, 13);
        graph.add_edge(n3, n4, 12);
        graph.add_edge(n3, n4, 125);
        graph.add_edge(n4, n5, 13);
        graph.add_edge(n5, n3, 13);
        debug_assert!(!is_strongly_connected(&graph));
        debug_assert_eq!(
            decompose_strongly_connected_components(&graph),
            vec![n0, n1, n1, n3, n3, n3]
        );
    }

    #[test]
    fn test_decompose_sccs_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        debug_assert!(is_strongly_connected(&graph));
        debug_assert_eq!(decompose_strongly_connected_components(&graph), vec![]);
    }

    #[test]
    fn test_decompose_sccs_nearly_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n1, 14);
        graph.add_edge(n1, n4, 15);
        graph.add_edge(n1, n4, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n0, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
        debug_assert_eq!(
            decompose_strongly_connected_components(&graph),
            vec![n0, n1, n1, n1, n1]
        );
    }

    #[test]
    fn test_decompose_sccs_nearly_scc_reverse() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n4, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n0, 14);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n1, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
        debug_assert_eq!(
            decompose_strongly_connected_components(&graph),
            vec![n0, n1, n1, n1, n1]
        );
    }

    /////////////////////////////////////////////////////////

    #[test]
    fn test_extract_subgraphs_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n0, 14);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n0, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        debug_assert_eq!(1, extracted.len());
        debug_assert!(is_strongly_connected(&extracted[0]));
        debug_assert_eq!(5, extracted[0].node_count());
        debug_assert_eq!(12, extracted[0].edge_count());
    }

    #[test]
    fn test_extract_subgraphs_one_wc_with_two_sccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 1);
        graph.add_edge(n0, n3, 18);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        debug_assert_eq!(2, extracted.len());
        for (i, graph) in extracted.iter().enumerate() {
            debug_assert!(
                is_strongly_connected(&extracted[i]),
                "Graph {} not strongly connected: {:?}",
                i,
                graph
            );
        }

        debug_assert_eq!(2, extracted[0].node_count());
        debug_assert_eq!(3, extracted[0].edge_count());
        debug_assert_eq!(3, extracted[1].node_count());
        debug_assert_eq!(5, extracted[1].edge_count());
    }

    #[test]
    fn test_extract_subgraphs_multiple_wccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        graph.add_edge(n0, n0, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n1, 13);
        graph.add_edge(n3, n4, 12);
        graph.add_edge(n3, n4, 125);
        graph.add_edge(n4, n5, 13);
        graph.add_edge(n5, n3, 13);
        debug_assert!(!is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        debug_assert_eq!(3, extracted.len());
        for (i, graph) in extracted.iter().enumerate() {
            debug_assert!(
                is_strongly_connected(&extracted[i]),
                "Graph {} not strongly connected: {:?}",
                i,
                graph
            );
        }

        debug_assert_eq!(1, extracted[0].node_count());
        debug_assert_eq!(1, extracted[0].edge_count());
        debug_assert_eq!(2, extracted[1].node_count());
        debug_assert_eq!(2, extracted[1].edge_count());
        debug_assert_eq!(3, extracted[2].node_count());
        debug_assert_eq!(4, extracted[2].edge_count());
    }

    #[test]
    fn test_extract_subgraphs_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        debug_assert!(is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        debug_assert!(extracted.is_empty());
    }

    #[test]
    fn test_extract_subgraphs_nearly_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n1, 14);
        graph.add_edge(n1, n4, 15);
        graph.add_edge(n1, n4, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n0, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        debug_assert_eq!(2, extracted.len());
        for (i, graph) in extracted.iter().enumerate() {
            debug_assert!(
                is_strongly_connected(&extracted[i]),
                "Graph {} not strongly connected: {:?}",
                i,
                graph
            );
        }

        debug_assert_eq!(1, extracted[0].node_count());
        debug_assert_eq!(0, extracted[0].edge_count());
        debug_assert_eq!(4, extracted[1].node_count());
        debug_assert_eq!(10, extracted[1].edge_count());
    }

    #[test]
    fn test_extract_subgraphs_nearly_scc_reverse() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        graph.add_edge(n4, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n4, n0, 14);
        graph.add_edge(n1, n0, 15);
        graph.add_edge(n1, n0, 155);
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n1, n4, 19);
        graph.add_edge(n2, n2, 20);
        debug_assert!(!is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        debug_assert_eq!(2, extracted.len());
        for (i, graph) in extracted.iter().enumerate() {
            debug_assert!(
                is_strongly_connected(&extracted[i]),
                "Graph {} not strongly connected: {:?}",
                i,
                graph
            );
        }

        debug_assert_eq!(1, extracted[0].node_count());
        debug_assert_eq!(0, extracted[0].edge_count());
        debug_assert_eq!(4, extracted[1].node_count());
        debug_assert_eq!(9, extracted[1].edge_count());
    }
}
