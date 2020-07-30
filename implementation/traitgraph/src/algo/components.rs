use crate::index::GraphIndex;
use crate::traversal::UndirectedBfs;
use crate::{MutableGraphContainer, StaticGraph};
use std::collections::LinkedList;

pub fn decompose_weakly_connected_components<Graph: Default + MutableGraphContainer + StaticGraph>(
    graph: &Graph,
) -> Vec<Graph>
where
    Graph::NodeData: Clone,
    Graph::EdgeData: Clone,
{
    let mut result = Vec::new();
    let mut nodes: LinkedList<_> = graph.node_indices().collect();
    // TODO this is not optimal. The Bfs recreates a vector of all nodes all the time.
    // Instead of doing that, the Bfs could reuse the order vector.
    // Then, an offset would need to be used for the subgraph indices.
    let mut visited = vec![false; graph.node_count()];

    while !nodes.is_empty() {
        let start = nodes.pop_front().unwrap();
        if visited[start.as_usize()] {
            continue;
        }

        let mut bfs: UndirectedBfs<Graph> = UndirectedBfs::new(graph, start);
        let mut subgraph = Graph::default();

        while let Some(node) = bfs.next(graph) {
            visited[node.as_usize()] = true;
            //println!("add_node: {:?}", node);
            subgraph.add_node(graph.node_data(node).clone());
            let subnode = bfs.order_of(node).unwrap();

            for out_neighbor in graph.out_neighbors(node) {
                let neighbor_id = out_neighbor.node_id;
                debug_assert!(
                    graph.contains_edge(node, neighbor_id),
                    "f: Edge missing: ({:?}, {:?})",
                    node,
                    neighbor_id
                );
                let edge_id = out_neighbor.edge_id;
                if let Some(subneighbor) = bfs.order_of(neighbor_id) {
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

                if let Some(subneighbor) = bfs.order_of(neighbor_id) {
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

#[cfg(test)]
mod tests {
    use crate::algo::components::decompose_weakly_connected_components;
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
}