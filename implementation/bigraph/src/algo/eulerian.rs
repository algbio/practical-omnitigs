use crate::interface::static_bigraph::StaticBigraph;
use crate::traitgraph::index::GraphIndex;
use bitvector::BitVector;
use std::collections::LinkedList;
use std::fmt::Write;
use traitgraph::walks::VecNodeWalk;

/*fn bitvector_to_string(bitvector: &BitVector) -> String {
    let mut result = String::new();
    for i in 0..bitvector.capacity() {
        if bitvector.contains(i) {
            write!(result, "1").unwrap()
        } else {
            write!(result, "0").unwrap()
        }
    }
    result
}*/

fn bitvector_to_index_string(bitvector: &BitVector) -> String {
    let mut result = String::new();
    for i in 0..bitvector.capacity() {
        if bitvector.contains(i) {
            write!(result, "{} ", i).unwrap()
        }
    }
    result
}

fn mark_edge_and_mirror<Graph: StaticBigraph>(
    bitvector: &mut BitVector,
    graph: &Graph,
    edge_index: Graph::EdgeIndex,
) {
    assert!(!bitvector.contains(edge_index.as_usize()));
    bitvector.insert(edge_index.as_usize());
    //println!("Marked edge {}", edge_index.as_usize());
    let mut inserted_mirror = false;
    let topological_mirror_edges = graph.topological_mirror_edges(edge_index);
    if topological_mirror_edges.contains(&edge_index) {
        inserted_mirror = true;
        //println!("Topological mirror edges contains self and is {:?}", topological_mirror_edges);
    } else {
        for mirror_edge in topological_mirror_edges {
            if bitvector.insert(mirror_edge.as_usize()) {
                //println!("Marked edge {}", mirror_edge.as_usize());
                inserted_mirror = true;
                break;
            }
        }
    }
    assert!(
        inserted_mirror,
        "bitvector {}\nforwards {:?}\nmirrors {:?}",
        bitvector_to_index_string(bitvector),
        graph
            .edges_between(
                graph.edge_endpoints(edge_index).from_node,
                graph.edge_endpoints(edge_index).to_node
            )
            .collect::<Vec<_>>(),
        graph.topological_mirror_edges(edge_index)
    );
}

/// Computes a Eulerian cycle in the graph.
pub fn compute_minimum_bidirected_eulerian_cycle_decomposition<Graph: StaticBigraph>(
    graph: &Graph,
) -> Vec<VecNodeWalk<Graph>> {
    let mut used_edges = BitVector::new(graph.edge_count());
    let mut cycles = Vec::new();

    for edge_index in graph.edge_indices() {
        if used_edges.contains(edge_index.as_usize()) {
            continue;
        }
        //println!("Starting new cycle with edge {}", edge_index.as_usize());

        let mut start_edge = Some(edge_index);
        let mut cycle = LinkedList::new();

        while let Some(start_edge_index) = start_edge {
            //println!("Start edge {}", start_edge_index.as_usize());
            mark_edge_and_mirror(&mut used_edges, graph, start_edge_index);
            let start_node = graph.edge_endpoints(start_edge_index).from_node;
            cycle.push_back(start_node);
            let mut current_node = graph.edge_endpoints(start_edge_index).to_node;

            while current_node != start_node {
                //println!("Expanding node {}", current_node.as_usize());
                cycle.push_back(current_node);
                let mut found_edge = false;
                for neighbor in graph.out_neighbors(current_node) {
                    if !used_edges.contains(neighbor.edge_id.as_usize()) {
                        mark_edge_and_mirror(&mut used_edges, graph, neighbor.edge_id);
                        found_edge = true;
                        current_node = graph.edge_endpoints(neighbor.edge_id).to_node;
                        break;
                    }
                }
                assert!(
                    found_edge,
                    "Found no continuation edge at node {}, d- {:?}, d+ {:?}\nbitvector {}",
                    current_node.as_usize(),
                    graph.in_neighbors(current_node).collect::<Vec<_>>(),
                    graph.out_neighbors(current_node).collect::<Vec<_>>(),
                    bitvector_to_index_string(&used_edges)
                );
            }

            //println!("Closed cycle, used_edges: {}", bitvector_to_string(&used_edges));

            // Find new start edge
            start_edge = None;
            for (cycle_index, &node_index) in cycle.iter().enumerate() {
                let mut found_neighbor = false;
                for neighbor in graph.out_neighbors(node_index) {
                    //println!("Found edge to continue current cycle {}", neighbor.edge_id.as_usize());
                    if !used_edges.contains(neighbor.edge_id.as_usize()) {
                        start_edge = Some(neighbor.edge_id);
                        found_neighbor = true;
                        break;
                    }
                }

                if found_neighbor {
                    //println!("Rotating at {}", cycle_index);
                    let mut rotator = cycle.split_off(cycle_index);
                    rotator.append(&mut cycle);
                    cycle = rotator;
                    //println!("Rotated cycle {:?}", cycle);
                    break;
                }
            }
        }

        let mut cycle_walk = Vec::new();
        cycle_walk.extend(cycle.iter());
        cycles.push(VecNodeWalk::new(cycle_walk));
    }

    cycles
}

/// Returns true if the graph contains a Eulerian bicycle.
pub fn decomposes_into_eulerian_bicycles<Graph: StaticBigraph>(graph: &Graph) -> bool {
    for node_index in graph.node_indices() {
        let incoming_self_loops = graph
            .in_neighbors(node_index)
            .filter(|n| graph.mirror_node(n.node_id).unwrap() == node_index)
            .count();
        let outgoing_self_loops = graph
            .out_neighbors(node_index)
            .filter(|n| graph.mirror_node(n.node_id).unwrap() == node_index)
            .count();

        if graph.in_degree(node_index) - incoming_self_loops
            != graph.out_degree(node_index) - outgoing_self_loops
        {
            return false;
        }
    }

    true
}

/// Compute a vector of nodes that has indegree != outdegree.
/// Bidirected self loops (edges from a node to its mirror) are handled by not counting them.
pub fn find_non_eulerian_binodes<Graph: StaticBigraph>(graph: &Graph) -> Vec<Graph::NodeIndex> {
    let mut result = Vec::new();
    for node_index in graph.node_indices() {
        let incoming_self_loops = graph
            .in_neighbors(node_index)
            .filter(|n| graph.mirror_node(n.node_id).unwrap() == node_index)
            .count();
        let outgoing_self_loops = graph
            .out_neighbors(node_index)
            .filter(|n| graph.mirror_node(n.node_id).unwrap() == node_index)
            .count();

        if graph.in_degree(node_index) - incoming_self_loops
            != graph.out_degree(node_index) - outgoing_self_loops
        {
            result.push(node_index);
        }
    }
    result
}

/// Compute a vector of tuples of nodes and outdegree - indegree that has indegree != outdegree.
/// Bidirected self loops (edges from a node to its mirror) are handled by not counting them.
pub fn find_non_eulerian_binodes_with_differences<Graph: StaticBigraph>(
    graph: &Graph,
) -> Vec<(Graph::NodeIndex, isize)> {
    let mut node_indices_and_differences = Vec::new();
    for node_index in graph.node_indices() {
        let incoming_self_loops = graph
            .in_neighbors(node_index)
            .filter(|n| graph.mirror_node(n.node_id).unwrap() == node_index)
            .count();
        let outgoing_self_loops = graph
            .out_neighbors(node_index)
            .filter(|n| graph.mirror_node(n.node_id).unwrap() == node_index)
            .count();

        let difference = (graph.out_degree(node_index) - outgoing_self_loops) as isize
            - (graph.in_degree(node_index) - incoming_self_loops) as isize;
        if difference != 0 {
            node_indices_and_differences.push((node_index, difference));
        }
    }
    node_indices_and_differences
}

#[cfg(test)]
mod tests {
    use crate::algo::eulerian::compute_minimum_bidirected_eulerian_cycle_decomposition;
    use crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
    use crate::interface::dynamic_bigraph::DynamicBigraph;
    use crate::interface::static_bigraph::StaticBigraphFromDigraph;
    use crate::traitgraph::interface::MutableGraphContainer;
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::walks::VecNodeWalk;

    #[test]
    fn test_bidirected_eulerian_cycle_triangle() {
        let mut bigraph = NodeBigraphWrapper::new(petgraph_impl::new());
        let n1 = bigraph.add_node(());
        let n2 = bigraph.add_node(());
        bigraph.set_mirror_nodes(n1, n2);
        let n3 = bigraph.add_node(());
        let n4 = bigraph.add_node(());
        bigraph.set_mirror_nodes(n3, n4);
        let n5 = bigraph.add_node(());
        let n6 = bigraph.add_node(());
        bigraph.set_mirror_nodes(n5, n6);

        bigraph.add_edge(n1, n3, 13);
        bigraph.add_edge(n4, n2, 42);
        bigraph.add_edge(n3, n5, 35);
        bigraph.add_edge(n6, n4, 64);
        bigraph.add_edge(n5, n1, 51);
        bigraph.add_edge(n2, n6, 26);

        let euler_cycles = compute_minimum_bidirected_eulerian_cycle_decomposition(&bigraph);
        assert_eq!(euler_cycles, vec![VecNodeWalk::new(vec![n1, n3, n5])]);
    }

    #[test]
    fn test_bidirected_eulerian_cycle_figure_eight() {
        let mut bigraph = NodeBigraphWrapper::new(petgraph_impl::new());
        let n1 = bigraph.add_node(());
        let n2 = bigraph.add_node(());
        bigraph.set_mirror_nodes(n1, n2);
        let n3 = bigraph.add_node(());
        let n4 = bigraph.add_node(());
        bigraph.set_mirror_nodes(n3, n4);
        let n5 = bigraph.add_node(());
        let n6 = bigraph.add_node(());
        bigraph.set_mirror_nodes(n5, n6);
        let n7 = bigraph.add_node(());
        let n8 = bigraph.add_node(());
        bigraph.set_mirror_nodes(n7, n8);
        let n9 = bigraph.add_node(());
        let n10 = bigraph.add_node(());
        bigraph.set_mirror_nodes(n9, n10);

        bigraph.add_edge(n1, n3, 13);
        bigraph.add_edge(n4, n2, 42);
        bigraph.add_edge(n3, n5, 35);
        bigraph.add_edge(n6, n4, 64);
        bigraph.add_edge(n5, n1, 51);
        bigraph.add_edge(n2, n6, 26);
        bigraph.add_edge(n1, n7, 13);
        bigraph.add_edge(n8, n2, 42);
        bigraph.add_edge(n7, n9, 35);
        bigraph.add_edge(n10, n8, 64);
        bigraph.add_edge(n9, n1, 51);
        bigraph.add_edge(n2, n10, 26);

        let euler_cycles = compute_minimum_bidirected_eulerian_cycle_decomposition(&bigraph);
        assert_eq!(
            euler_cycles,
            vec![VecNodeWalk::new(vec![n1, n3, n5, n1, n7, n9])]
        );
    }
}
