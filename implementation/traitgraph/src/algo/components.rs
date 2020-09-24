use super::traversal::{
    AllowedNodesForbiddenSubgraph, PostOrderForwardDfs, PreOrderBackwardBfs, PreOrderForwardBfs,
    PreOrderUndirectedBfs,
};
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
    // Using a vector might be faster here. We anyways only pop from one side.
    let mut nodes: LinkedList<_> = graph.node_indices().collect();
    // TODO this is not optimal. The bfs recreates a vector of all nodes all the time.
    // Instead of doing that, the bfs could reuse the order vector.
    // Then, an offset would need to be used for the subgraph indices.
    let mut visited = vec![false; graph.node_count()];

    while !nodes.is_empty() {
        let start = nodes.pop_front().unwrap();
        if visited[start.as_usize()] {
            continue;
        }

        let mut bfs = PreOrderUndirectedBfs::new(graph, start);
        let mut subgraph = Graph::default();

        while let Some(node) = bfs.next() {
            let node = if let NodeOrEdge::Node(node) = node {
                node
            } else {
                continue;
            };
            visited[node.as_usize()] = true;
            //println!("add_node: {:?}", node);
            subgraph.add_node(graph.node_data(node).clone());
            let subnode = bfs.rank_of(node).unwrap();

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
    // TODO this is not optimal. The dfs recreates a vector of all nodes all the time.
    // Instead of doing that, the dfs could reuse the order vector.
    // Then, an offset would need to be used for the subgraph indices.
    let mut visited = vec![false; graph.node_count()];

    for node in graph.node_indices() {
        if !visited[node.as_usize()] {
            let mut dfs = PostOrderForwardDfs::new(graph, node);

            while let Some(node) = dfs.next(graph) {
                visited[node.as_usize()] = true;
                nodes.push_front(node);
            }
        }
    }

    //println!("nodes: {:?}", nodes);

    for root_node in nodes {
        if visited[root_node.as_usize()] {
            //println!("Reverse processing {:?}", root_node);
            let mut bfs = PreOrderBackwardBfs::new(graph, root_node);

            while let Some(node) =
                bfs.next_with_forbidden_subgraph(&AllowedNodesForbiddenSubgraph::new(&visited))
            {
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
        let e0 = graph.add_edge(n0, n1, 10);
        let e1 = graph.add_edge(n1, n2, 11);
        let e15 = graph.add_edge(n1, n2, 115);
        let e2 = graph.add_edge(n2, n3, 12);
        let e3 = graph.add_edge(n3, n4, 13);
        let e4 = graph.add_edge(n4, n0, 14);
        let e45 = graph.add_edge(n4, n0, 145);
        let e5 = graph.add_edge(n1, n0, 15);
        let e6 = graph.add_edge(n2, n1, 16);
        let e7 = graph.add_edge(n3, n2, 17);
        let e8 = graph.add_edge(n4, n3, 18);
        let e9 = graph.add_edge(n0, n4, 19);
        let e10 = graph.add_edge(n2, n2, 20);
        let result = decompose_weakly_connected_components(&graph);
        assert_eq!(result.len(), 1);
        let result = result.first().unwrap();
        assert_eq!(result.node_count(), graph.node_count());
        assert_eq!(result.edge_count(), graph.edge_count());

        assert_eq!(result.node_data(n0), &0);
        assert_eq!(result.node_data(n1), &4);
        assert_eq!(result.node_data(n2), &1);
        assert_eq!(result.node_data(n3), &3);
        assert_eq!(result.node_data(n4), &2);

        assert_eq!(result.edge_data(e0), &145);
        assert_eq!(result.edge_data(e1), &14);
        assert_eq!(result.edge_data(e15), &19);
        assert_eq!(result.edge_data(e2), &15);
        assert_eq!(result.edge_data(e3), &10);
        assert_eq!(result.edge_data(e4), &13);
        assert_eq!(result.edge_data(e45), &18);
        assert_eq!(result.edge_data(e5), &20);
        assert_eq!(result.edge_data(e6), &16);
        assert_eq!(result.edge_data(e7), &12);
        assert_eq!(result.edge_data(e8), &17);
        assert_eq!(result.edge_data(e9), &115);
        assert_eq!(result.edge_data(e10), &11);
    }

    #[test]
    fn test_decompose_weakly_connected_components_one_wc_with_two_sccs() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let e0 = graph.add_edge(n0, n1, 10);
        let e1 = graph.add_edge(n1, n2, 11);
        let e2 = graph.add_edge(n2, n3, 12);
        let e3 = graph.add_edge(n3, n4, 13);
        let e4 = graph.add_edge(n1, n0, 15);
        let e5 = graph.add_edge(n3, n2, 17);
        let e6 = graph.add_edge(n4, n3, 18);
        let e7 = graph.add_edge(n2, n2, 20);
        let e8 = graph.add_edge(n2, n2, 21);
        let result = decompose_weakly_connected_components(&graph);
        assert_eq!(result.len(), 1);
        let result = result.first().unwrap();
        assert_eq!(result.node_count(), graph.node_count());
        assert_eq!(result.edge_count(), graph.edge_count());

        assert_eq!(result.node_data(n0), &0);
        assert_eq!(result.node_data(n1), &1);
        assert_eq!(result.node_data(n2), &2);
        assert_eq!(result.node_data(n3), &3);
        assert_eq!(result.node_data(n4), &4);

        assert_eq!(result.edge_data(e0), &15);
        assert_eq!(result.edge_data(e1), &10);
        assert_eq!(result.edge_data(e2), &21);
        assert_eq!(result.edge_data(e3), &20);
        assert_eq!(result.edge_data(e4), &11);
        assert_eq!(result.edge_data(e5), &17);
        assert_eq!(result.edge_data(e6), &12);
        assert_eq!(result.edge_data(e7), &18);
        assert_eq!(result.edge_data(e8), &13);
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
        assert_eq!(result.len(), 3);
        let first = result[0].clone();
        let second = result[1].clone();
        let third = result[2].clone();
        assert_eq!(first.node_count(), 1);
        assert_eq!(first.edge_count(), 1);
        assert_eq!(second.node_count(), 2);
        assert_eq!(second.edge_count(), 2);
        assert_eq!(third.node_count(), 3);
        assert_eq!(third.edge_count(), 2);

        assert_eq!(first.node_data(0.into()), &0);
        assert_eq!(second.node_data(0.into()), &1);
        assert_eq!(second.node_data(1.into()), &2);
        assert_eq!(third.node_data(0.into()), &3);
        assert_eq!(third.node_data(1.into()), &4);
        assert_eq!(third.node_data(2.into()), &5);

        assert_eq!(first.edge_data(0.into()), &10);
        assert_eq!(second.edge_data(0.into()), &115);
        assert_eq!(second.edge_data(1.into()), &11);
        assert_eq!(third.edge_data(0.into()), &12);
        assert_eq!(third.edge_data(1.into()), &13);
    }

    #[test]
    fn test_decompose_weakly_connected_components_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        let result = decompose_weakly_connected_components(&graph);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_decompose_weakly_connected_components_nearly_scc() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let e0 = graph.add_edge(n0, n1, 10);
        let e05 = graph.add_edge(n0, n1, 105);
        let e1 = graph.add_edge(n1, n2, 11);
        let e2 = graph.add_edge(n2, n3, 12);
        let e3 = graph.add_edge(n3, n4, 13);
        let e4 = graph.add_edge(n4, n1, 14);
        let e5 = graph.add_edge(n1, n4, 15);
        let e6 = graph.add_edge(n2, n1, 16);
        let e7 = graph.add_edge(n3, n2, 17);
        let e8 = graph.add_edge(n4, n3, 18);
        let e9 = graph.add_edge(n0, n4, 19);
        let e10 = graph.add_edge(n2, n2, 20);
        let result = decompose_weakly_connected_components(&graph);
        assert_eq!(result.len(), 1);
        let result = result.first().unwrap();
        assert_eq!(result.node_count(), graph.node_count());
        assert_eq!(result.edge_count(), graph.edge_count());

        assert_eq!(result.node_data(n0), &0);
        assert_eq!(result.node_data(n1), &4);
        assert_eq!(result.node_data(n2), &1);
        assert_eq!(result.node_data(n3), &3);
        assert_eq!(result.node_data(n4), &2);

        assert_eq!(result.edge_data(e0), &19);
        assert_eq!(result.edge_data(e05), &15);
        assert_eq!(result.edge_data(e1), &14);
        assert_eq!(result.edge_data(e2), &105);
        assert_eq!(result.edge_data(e3), &10);
        assert_eq!(result.edge_data(e4), &13);
        assert_eq!(result.edge_data(e5), &18);
        assert_eq!(result.edge_data(e6), &20);
        assert_eq!(result.edge_data(e7), &16);
        assert_eq!(result.edge_data(e8), &12);
        assert_eq!(result.edge_data(e9), &17);
        assert_eq!(result.edge_data(e10), &11);
    }

    #[test]
    fn test_decompose_weakly_connected_components_nearly_scc_reverse() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let e0 = graph.add_edge(n4, n1, 10);
        let e1 = graph.add_edge(n1, n2, 11);
        let e2 = graph.add_edge(n2, n3, 12);
        let e3 = graph.add_edge(n3, n4, 13);
        let e4 = graph.add_edge(n4, n0, 14);
        let e5 = graph.add_edge(n1, n0, 15);
        let e55 = graph.add_edge(n1, n0, 155);
        let e6 = graph.add_edge(n2, n1, 16);
        let e7 = graph.add_edge(n3, n2, 17);
        let e8 = graph.add_edge(n4, n3, 18);
        let e9 = graph.add_edge(n1, n4, 19);
        let e10 = graph.add_edge(n2, n2, 20);
        let result = decompose_weakly_connected_components(&graph);
        assert_eq!(result.len(), 1);
        let result = result.first().unwrap();
        assert_eq!(result.node_count(), graph.node_count());
        assert_eq!(result.edge_count(), graph.edge_count());

        assert_eq!(result.node_data(n0), &0);
        assert_eq!(result.node_data(n1), &1);
        assert_eq!(result.node_data(n2), &4);
        assert_eq!(result.node_data(n3), &2);
        assert_eq!(result.node_data(n4), &3);

        assert_eq!(result.edge_data(e0), &155);
        assert_eq!(result.edge_data(e1), &15);
        assert_eq!(result.edge_data(e2), &14);
        assert_eq!(result.edge_data(e3), &10);
        assert_eq!(result.edge_data(e4), &19);
        assert_eq!(result.edge_data(e5), &20);
        assert_eq!(result.edge_data(e55), &16);
        assert_eq!(result.edge_data(e6), &11);
        assert_eq!(result.edge_data(e7), &17);
        assert_eq!(result.edge_data(e8), &13);
        assert_eq!(result.edge_data(e9), &18);
        assert_eq!(result.edge_data(e10), &12);
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
        assert!(is_strongly_connected(&graph));
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
        assert!(!is_strongly_connected(&graph));
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
        assert!(!is_strongly_connected(&graph));
    }

    #[test]
    fn test_scc_check_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        assert!(is_strongly_connected(&graph));
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
        assert!(!is_strongly_connected(&graph));
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
        assert!(!is_strongly_connected(&graph));
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
        assert!(is_strongly_connected(&graph));
        assert_eq!(decompose_strongly_connected_components(&graph), vec![n0; 5]);
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
        assert!(!is_strongly_connected(&graph));
        assert_eq!(
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
        assert!(!is_strongly_connected(&graph));
        assert_eq!(
            decompose_strongly_connected_components(&graph),
            vec![n0, n1, n1, n3, n3, n3]
        );
    }

    #[test]
    fn test_decompose_sccs_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        assert!(is_strongly_connected(&graph));
        assert_eq!(decompose_strongly_connected_components(&graph), vec![]);
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
        assert!(!is_strongly_connected(&graph));
        assert_eq!(
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
        assert!(!is_strongly_connected(&graph));
        assert_eq!(
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
        assert!(is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        assert_eq!(1, extracted.len());
        assert!(is_strongly_connected(&extracted[0]));
        assert_eq!(5, extracted[0].node_count());
        assert_eq!(12, extracted[0].edge_count());
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
        assert!(!is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        assert_eq!(2, extracted.len());
        for (i, graph) in extracted.iter().enumerate() {
            assert!(
                is_strongly_connected(&extracted[i]),
                "Graph {} not strongly connected: {:?}",
                i,
                graph
            );
        }

        assert_eq!(2, extracted[0].node_count());
        assert_eq!(3, extracted[0].edge_count());
        assert_eq!(3, extracted[1].node_count());
        assert_eq!(5, extracted[1].edge_count());
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
        assert!(!is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        assert_eq!(3, extracted.len());
        for (i, graph) in extracted.iter().enumerate() {
            assert!(
                is_strongly_connected(&extracted[i]),
                "Graph {} not strongly connected: {:?}",
                i,
                graph
            );
        }

        assert_eq!(1, extracted[0].node_count());
        assert_eq!(1, extracted[0].edge_count());
        assert_eq!(2, extracted[1].node_count());
        assert_eq!(2, extracted[1].edge_count());
        assert_eq!(3, extracted[2].node_count());
        assert_eq!(4, extracted[2].edge_count());
    }

    #[test]
    fn test_extract_subgraphs_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        assert!(is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        assert!(extracted.is_empty());
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
        assert!(!is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        assert_eq!(2, extracted.len());
        for (i, graph) in extracted.iter().enumerate() {
            assert!(
                is_strongly_connected(&extracted[i]),
                "Graph {} not strongly connected: {:?}",
                i,
                graph
            );
        }

        assert_eq!(1, extracted[0].node_count());
        assert_eq!(0, extracted[0].edge_count());
        assert_eq!(4, extracted[1].node_count());
        assert_eq!(10, extracted[1].edge_count());
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
        assert!(!is_strongly_connected(&graph));
        let extracted = extract_subgraphs_from_node_mapping(
            &graph,
            &decompose_strongly_connected_components(&graph),
        );
        assert_eq!(2, extracted.len());
        for (i, graph) in extracted.iter().enumerate() {
            assert!(
                is_strongly_connected(&extracted[i]),
                "Graph {} not strongly connected: {:?}",
                i,
                graph
            );
        }

        assert_eq!(1, extracted[0].node_count());
        assert_eq!(0, extracted[0].edge_count());
        assert_eq!(4, extracted[1].node_count());
        assert_eq!(9, extracted[1].edge_count());
    }
}
