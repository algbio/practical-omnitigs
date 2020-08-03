use crate::index::GraphIndex;
use crate::traversal::{PreOrderBackwardBfs, PreOrderForwardBfs, PreOrderUndirectedBfs};
use crate::{MutableGraphContainer, StaticGraph};
use std::collections::LinkedList;

/// Returns the weakly connected components of a graph.
///
/// If the graph is empty, no WCCs are returned.
/// Otherwise, the WCCs are cloned into new graphs, without preserving the node or edge indices.
pub fn decompose_weakly_connected_components<Graph: Default + MutableGraphContainer + StaticGraph>(
    graph: &Graph,
) -> Vec<Graph>
where
    Graph::NodeData: Clone,
    Graph::EdgeData: Clone,
{
    let mut result = Vec::new();
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

        while let Some(node) = bfs.next(graph) {
            visited[node.as_usize()] = true;
            //println!("add_node: {:?}", node);
            subgraph.add_node(graph.node_data(node).clone());
            let subnode = bfs.rank_of(node).unwrap();

            for out_neighbor in graph.out_neighbors(node) {
                let neighbor_id = out_neighbor.node_id;
                debug_assert!(
                    graph.contains_edge(node, neighbor_id),
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
                    graph.contains_edge(neighbor_id, node),
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

    let mut traversal = PreOrderForwardBfs::new(graph, graph.node_indices().next().unwrap());
    let mut traversal_node_count = 0;

    while traversal.next(graph).is_some() {
        traversal_node_count += 1;
    }

    if traversal_node_count != graph.node_count() {
        debug_assert!(traversal_node_count < graph.node_count());
        return false;
    }

    let mut traversal = PreOrderBackwardBfs::new(graph, graph.node_indices().next().unwrap());
    let mut traversal_node_count = 0;

    while traversal.next(graph).is_some() {
        traversal_node_count += 1;
    }

    if traversal_node_count != graph.node_count() {
        debug_assert!(traversal_node_count < graph.node_count());
        return false;
    }

    true
}

#[cfg(test)]
mod tests {
    use crate::algo::components::{decompose_weakly_connected_components, is_strongly_connected};
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
        let e2 = graph.add_edge(n2, n3, 12);
        let e3 = graph.add_edge(n3, n4, 13);
        let e4 = graph.add_edge(n4, n0, 14);
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

        assert_eq!(result.edge_data(e0), &14);
        assert_eq!(result.edge_data(e1), &19);
        assert_eq!(result.edge_data(e2), &15);
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
        assert_eq!(result.edge_data(e2), &20);
        assert_eq!(result.edge_data(e3), &11);
        assert_eq!(result.edge_data(e4), &17);
        assert_eq!(result.edge_data(e5), &12);
        assert_eq!(result.edge_data(e6), &18);
        assert_eq!(result.edge_data(e7), &13);
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
        assert_eq!(second.edge_count(), 1);
        assert_eq!(third.node_count(), 3);
        assert_eq!(third.edge_count(), 2);

        assert_eq!(first.node_data(0.into()), &0);
        assert_eq!(second.node_data(0.into()), &1);
        assert_eq!(second.node_data(1.into()), &2);
        assert_eq!(third.node_data(0.into()), &3);
        assert_eq!(third.node_data(1.into()), &4);
        assert_eq!(third.node_data(2.into()), &5);

        assert_eq!(first.edge_data(0.into()), &10);
        assert_eq!(second.edge_data(0.into()), &11);
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
        assert_eq!(result.edge_data(e1), &15);
        assert_eq!(result.edge_data(e2), &14);
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

        assert_eq!(result.edge_data(e0), &15);
        assert_eq!(result.edge_data(e1), &14);
        assert_eq!(result.edge_data(e2), &10);
        assert_eq!(result.edge_data(e3), &19);
        assert_eq!(result.edge_data(e4), &20);
        assert_eq!(result.edge_data(e5), &16);
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
        graph.add_edge(n3, n4, 12);
        graph.add_edge(n4, n5, 13);
        assert!(!is_strongly_connected(&graph));
    }

    #[test]
    fn test_scc_check_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();
        assert!(is_strongly_connected(&graph));
    }

    #[test]
    fn ttest_scc_check_nearly_scc() {
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
        graph.add_edge(n2, n1, 16);
        graph.add_edge(n3, n2, 17);
        graph.add_edge(n4, n3, 18);
        graph.add_edge(n1, n4, 19);
        graph.add_edge(n2, n2, 20);
        assert!(!is_strongly_connected(&graph));
    }
}
