use crate::index::GraphIndex;
use crate::interface::StaticGraph;

/// Returns the amount of unitigs in the graph that are longer than two nodes.
/// Those are called uncompacted, because they can be trivially compacted by contracting the inner nodes.
pub fn count_uncompacted_unitigs<Graph: StaticGraph>(graph: &Graph) -> usize {
    let mut used_nodes = vec![false; graph.node_count()];
    let mut uncompacted_unitig_count = 0;

    for node_index in graph.node_indices() {
        if used_nodes[node_index.as_usize()] {
            continue;
        } else {
            used_nodes[node_index.as_usize()] = true;
        }

        if graph.out_neighbors(node_index).into_iter().count() == 1
            && graph.in_neighbors(node_index).into_iter().count() == 1
        {
            uncompacted_unitig_count += 1;

            let mut start_index = node_index;
            let mut end_index = node_index;

            while graph.out_neighbors(start_index).into_iter().count() <= 1
                && graph.in_neighbors(start_index).into_iter().count() == 1
            {
                used_nodes[start_index.as_usize()] = true;
                start_index = graph
                    .in_neighbors(start_index)
                    .into_iter()
                    .next()
                    .unwrap()
                    .node_id;
            }
            while graph.out_neighbors(end_index).into_iter().count() == 1
                && graph.in_neighbors(end_index).into_iter().count() <= 1
            {
                used_nodes[end_index.as_usize()] = true;
                end_index = graph
                    .out_neighbors(end_index)
                    .into_iter()
                    .next()
                    .unwrap()
                    .node_id;
            }

            // println!("Found uncompacted unitig from {:?} to {:?}", start_index, end_index);
        }
    }

    uncompacted_unitig_count
}

#[cfg(test)]
mod tests {
    use crate::petgraph_impl;
    use crate::unitigs::count_uncompacted_unitigs;
    use crate::MutableGraphContainer;

    #[test]
    fn test_count_uncompacted_unitigs_simple() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let n6 = graph.add_node(6);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n0, n2, 11);
        graph.add_edge(n1, n2, 12);
        graph.add_edge(n2, n3, 13);
        graph.add_edge(n2, n4, 14);
        graph.add_edge(n3, n0, 15);
        graph.add_edge(n3, n4, 16);
        graph.add_edge(n4, n5, 17);
        graph.add_edge(n5, n6, 18);
        graph.add_edge(n6, n0, 19);
        assert_eq!(count_uncompacted_unitigs(&graph), 2);
    }
}
