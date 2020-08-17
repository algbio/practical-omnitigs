use crate::index::GraphIndex;
use crate::interface::StaticGraph;

#[derive(Debug, Eq, PartialEq, Default)]
pub struct UncompactedNodeUnitigs {
    /// Uncompacted unitigs of length two
    pub len_2: usize,
    /// Uncompacted unitigs of length three
    pub len_3: usize,
    /// Uncompacted unitigs of length four or more
    pub len_4_more: usize,
}

impl UncompactedNodeUnitigs {
    pub fn total(&self) -> usize {
        self.len_2 + self.len_3 + self.len_4_more
    }
}

/// Returns the amount of node unitigs in the graph that can be compacted.
/// A unitig can be compacted if
///  * it contains more than one inner node, or
///  * it contains one inner node and its first node has outdegree one, or
///  * it contains one inner node and its last node has indegree one, or
///  * it contains two nodes (no inner node) and its first node has outdegree one and its second node has indegree one.
pub fn count_uncompacted_node_unitigs<Graph: StaticGraph>(graph: &Graph) -> UncompactedNodeUnitigs {
    let mut used_edges = vec![false; graph.edge_count()];
    let mut uncompacted_unitig_count = UncompactedNodeUnitigs::default();

    for first_index in graph.node_indices() {
        for neighbor in graph.out_neighbors(first_index) {
            let edge_index = neighbor.edge_id;
            let second_index = neighbor.node_id;
            if used_edges[edge_index.as_usize()] {
                continue;
            } else {
                used_edges[edge_index.as_usize()] = true;
            }

            let mut start_index = first_index;
            let mut end_index = second_index;
            let mut length = 2usize;

            while graph.out_degree(start_index) == 1 && graph.in_degree(start_index) == 1 {
                let neighbor = graph.in_neighbors(start_index).into_iter().next().unwrap();

                start_index = neighbor.node_id;
                used_edges[neighbor.edge_id.as_usize()] = true;
                length += 1;
            }
            while graph.out_degree(end_index) == 1 && graph.in_degree(end_index) == 1 {
                let neighbor = graph.out_neighbors(end_index).into_iter().next().unwrap();

                end_index = neighbor.node_id;
                used_edges[neighbor.edge_id.as_usize()] = true;
                length += 1;
            }

            let start_out_degree = graph.out_degree(start_index);
            let end_in_degree = graph.in_degree(end_index);

            if length >= 4 {
                uncompacted_unitig_count.len_4_more += 1;
            //println!("Found uncompacted unitig from {:?} to {:?} (len: {})", start_index, end_index, length);
            } else if length >= 3 && (start_out_degree == 1 || end_in_degree == 1) {
                uncompacted_unitig_count.len_3 += 1;
            } else if length >= 2 && start_out_degree == 1 && end_in_degree == 1 {
                uncompacted_unitig_count.len_2 += 1;
            } else {
                //println!("Found compacted unitig from {:?} to {:?}", start_index, end_index);
            }
        }
    }

    uncompacted_unitig_count
}

#[cfg(test)]
mod tests {
    use super::count_uncompacted_node_unitigs;
    use super::UncompactedNodeUnitigs;
    use crate::implementation::petgraph_impl;
    use crate::interface::MutableGraphContainer;

    #[test]
    fn test_count_uncompacted_node_unitigs_four() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let n6 = graph.add_node(6);
        let n7 = graph.add_node(7);
        let n8 = graph.add_node(8);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n3, n5, 14);
        graph.add_edge(n4, n8, 15);
        graph.add_edge(n5, n8, 16);
        graph.add_edge(n8, n6, 17);
        graph.add_edge(n8, n7, 18);
        graph.add_edge(n6, n0, 19);
        graph.add_edge(n7, n0, 20);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 0,
                len_4_more: 1
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n8, n3, 21);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 0,
                len_4_more: 1
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n0, n8, 22);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 0,
                len_4_more: 1
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n8, n3, 21);
        graph2.add_edge(n0, n8, 22);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 0,
                len_4_more: 1
            }
        );

        drop(graph); // Against linter errors for last clone.
    }

    #[test]
    fn test_count_uncompacted_node_unitigs_three() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let n6 = graph.add_node(6);
        let n7 = graph.add_node(7);
        let n8 = graph.add_node(8);
        graph.add_edge(n0, n2, 10);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n3, n5, 14);
        graph.add_edge(n4, n8, 15);
        graph.add_edge(n5, n8, 16);
        graph.add_edge(n8, n6, 17);
        graph.add_edge(n8, n7, 18);
        graph.add_edge(n6, n0, 19);
        graph.add_edge(n7, n0, 20);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 1,
                len_4_more: 0
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n8, n3, 21);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 1,
                len_4_more: 0
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n0, n8, 22);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 1,
                len_4_more: 0
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n8, n3, 21);
        graph2.add_edge(n0, n8, 22);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 0,
                len_4_more: 0
            }
        );

        drop(graph); // Against linter errors for last clone.
    }

    #[test]
    fn test_count_uncompacted_node_unitigs_two() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let n6 = graph.add_node(6);
        let n7 = graph.add_node(7);
        let n8 = graph.add_node(8);
        graph.add_edge(n0, n3, 10);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n3, n5, 14);
        graph.add_edge(n4, n8, 15);
        graph.add_edge(n5, n8, 16);
        graph.add_edge(n8, n6, 17);
        graph.add_edge(n8, n7, 18);
        graph.add_edge(n6, n0, 19);
        graph.add_edge(n7, n0, 20);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph),
            UncompactedNodeUnitigs {
                len_2: 1,
                len_3: 0,
                len_4_more: 0
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n8, n3, 21);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 0,
                len_4_more: 0
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n0, n8, 22);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 0,
                len_4_more: 0
            }
        );

        let mut graph2 = graph.clone();
        graph2.add_edge(n8, n3, 21);
        graph2.add_edge(n0, n8, 22);
        assert_eq!(
            count_uncompacted_node_unitigs(&graph2),
            UncompactedNodeUnitigs {
                len_2: 0,
                len_3: 0,
                len_4_more: 0
            }
        );

        drop(graph); // Against linter errors for last clone.
    }
}
