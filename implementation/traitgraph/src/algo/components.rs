use crate::traversal::Bfs;
use crate::{ImmutableGraphContainer, MutableGraphContainer, NavigableGraph};
use num_traits::{NumCast, PrimInt};
use std::collections::LinkedList;

pub fn decompose_weakly_connected_components<
    'a,
    NodeData: Clone,
    EdgeData: Clone,
    IndexType: PrimInt,
    Graph: Default
        + MutableGraphContainer<NodeData, EdgeData, IndexType>
        + ImmutableGraphContainer<NodeData, EdgeData, IndexType>
        + NavigableGraph<'a, NodeData, EdgeData, IndexType>,
>(
    graph: &'a Graph,
) -> Vec<Graph> {
    let mut result = Vec::new();
    let mut nodes: LinkedList<_> = graph.node_indices().collect();
    // TODO this is not optimal. The Bfs recreates a vector of all nodes all the time.
    // Instead of doing that, the Bfs could reuse the order vector.
    // Then, an offset would need to be used for the subgraph indices.
    let mut visited = vec![false; graph.node_count()];

    while !nodes.is_empty() {
        let start = nodes.pop_front().unwrap();
        if visited[<usize as NumCast>::from(start).unwrap()] {
            continue;
        }

        let mut bfs = Bfs::new(graph, start);
        let mut subgraph = Graph::default();
        subgraph.add_node(graph.node_data(start).unwrap().clone());

        while let Some(node) = bfs.next(graph) {
            visited[<usize as NumCast>::from(node).unwrap()] = true;
            subgraph.add_node(graph.node_data(node).unwrap().clone());
            let subnode = bfs.order_of(node).unwrap().into();

            for out_neighbor in graph.out_neighbors(node).unwrap() {
                let neighbor_id = out_neighbor.node_id;
                let edge_id = out_neighbor.edge_id;
                if let Some(subneighbor) = bfs.order_of(neighbor_id).map(Into::into) {
                    let edge_data = graph.edge_data(edge_id).unwrap().clone();
                    subgraph.add_edge(subnode, subneighbor, edge_data);
                }
            }
            for in_neighbor in graph.in_neighbors(node).unwrap() {
                let neighbor_id = in_neighbor.node_id;
                let edge_id = in_neighbor.edge_id;
                if let Some(subneighbor) = bfs.order_of(neighbor_id).map(Into::into) {
                    let edge_data = graph.edge_data(edge_id).unwrap().clone();
                    subgraph.add_edge(subneighbor, subnode, edge_data);
                }
            }
        }

        result.push(subgraph);
    }

    result
}
