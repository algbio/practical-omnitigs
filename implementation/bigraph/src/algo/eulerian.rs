use crate::algo::{bitvector_to_index_string, mark_edge_and_mirror};
use crate::interface::static_bigraph::StaticEdgeCentricBigraph;
use crate::interface::BidirectedData;
use crate::traitgraph::index::GraphIndex;
use bitvector::BitVector;
use std::collections::{HashSet, LinkedList};
use traitgraph::walks::VecEdgeWalk;

/// Computes a Eulerian cycle in the graph.
pub fn compute_minimum_bidirected_eulerian_cycle_decomposition<
    EdgeData: BidirectedData + Eq,
    Graph: StaticEdgeCentricBigraph<EdgeData = EdgeData>,
>(
    graph: &Graph,
) -> Vec<VecEdgeWalk<Graph>> {
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
            cycle.push_back(start_edge_index);
            let mut current_node = graph.edge_endpoints(start_edge_index).to_node;

            let mut has_neighbor = true;
            while has_neighbor {
                //println!("Expanding node {}", current_node.as_usize());
                has_neighbor = false;
                for neighbor in graph.out_neighbors(current_node) {
                    if !used_edges.contains(neighbor.edge_id.as_usize()) {
                        cycle.push_back(neighbor.edge_id);
                        mark_edge_and_mirror(&mut used_edges, graph, neighbor.edge_id);
                        has_neighbor = true;
                        current_node = neighbor.node_id;
                        break;
                    }
                }
                debug_assert!(
                    has_neighbor || current_node == start_node,
                    "Found no continuation edge at node {}, d- {:?}, d+ {:?}, diff {}, self-mirror {}, #d- {:?}, #d+ {:?}\nbitvector {}",
                    current_node.as_usize(),
                    graph.in_neighbors(current_node).collect::<Vec<_>>(),
                    graph.out_neighbors(current_node).collect::<Vec<_>>(),
                    compute_eulerian_superfluous_out_biedges(graph, current_node),
                    graph.is_self_mirror_node(current_node),
                    graph.in_bidegree(current_node),
                    graph.out_bidegree(current_node),
                    bitvector_to_index_string(&used_edges),
                );
            }

            //println!("Closed cycle, used_edges: {}", bitvector_to_index_string(&used_edges));

            // Find new start edge
            start_edge = None;
            for (cycle_index, &edge_index) in cycle.iter().enumerate() {
                let mut found_neighbor = false;
                let from_node = graph.edge_endpoints(edge_index).from_node;
                for neighbor in graph.out_neighbors(from_node) {
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
        cycles.push(cycle_walk);
    }

    cycles
}

/// Returns true if the graph contains a Eulerian bicycle.
pub fn decomposes_into_eulerian_bicycles<
    EdgeData: BidirectedData + Eq,
    Graph: StaticEdgeCentricBigraph<EdgeData = EdgeData>,
>(
    graph: &Graph,
) -> bool {
    for node_index in graph.node_indices() {
        if compute_eulerian_superfluous_out_biedges(graph, node_index) != 0 {
            return false;
        }
    }

    true
}

/// Compute a vector of nodes that has indegree != outdegree.
/// Self-mirror nodes missing a biedge are returned as well (even though they always have indegree == outdegree).
pub fn find_non_eulerian_binodes<
    EdgeData: BidirectedData + Eq,
    Graph: StaticEdgeCentricBigraph<EdgeData = EdgeData>,
>(
    graph: &Graph,
) -> Vec<Graph::NodeIndex> {
    find_non_eulerian_binodes_with_differences(graph)
        .iter()
        .map(|&(n, _)| n)
        .collect()
}

/// Computes the number of outgoing edges that need to be added to make the binode Eulerian.
/// A negative number indicates that incoming edges need to be added.
pub fn compute_eulerian_superfluous_out_biedges<
    EdgeData: BidirectedData + Eq,
    Graph: StaticEdgeCentricBigraph<EdgeData = EdgeData>,
>(
    graph: &Graph,
    node_index: Graph::NodeIndex,
) -> isize {
    let mirror_node = graph.mirror_node(node_index).unwrap();
    if mirror_node == node_index {
        (graph.out_degree(node_index) % 2) as isize
    } else {
        let mut out_neighbor_count = 0;
        let mut out_inversion_count = 0;
        let mut out_inversion_mirrors = HashSet::new();
        for out_neighbor in graph.out_neighbors(node_index) {
            if out_neighbor.node_id == mirror_node {
                if out_inversion_mirrors.insert(out_neighbor.edge_id) {
                    out_inversion_count += 1;
                    let mirror_edge = graph
                        .mirror_edge_edge_centric(out_neighbor.edge_id)
                        .unwrap();
                    let inserted = out_inversion_mirrors.insert(mirror_edge);
                    debug_assert!(inserted);
                }
            } else {
                out_neighbor_count += 1;
            }
        }
        let mut in_neighbor_count = 0;
        let mut in_inversion_count = 0;
        let mut in_inversion_mirrors = HashSet::new();
        for in_neighbor in graph.in_neighbors(node_index) {
            if in_neighbor.node_id == mirror_node {
                if in_inversion_mirrors.insert(in_neighbor.edge_id) {
                    in_inversion_count += 1;
                    let mirror_edge = graph.mirror_edge_edge_centric(in_neighbor.edge_id).unwrap();
                    let inserted = in_inversion_mirrors.insert(mirror_edge);
                    debug_assert!(inserted);
                }
            } else {
                in_neighbor_count += 1;
            }
        }
        let inversion_diff = out_inversion_count - in_inversion_count;

        out_neighbor_count += inversion_diff;
        in_neighbor_count -= inversion_diff;
        out_neighbor_count - in_neighbor_count
    }
}

/// Compute a vector of tuples of nodes and outdegree - indegree that has indegree != outdegree.
/// Self-mirror nodes missing a biedge are returned with a 0 in the second entry of the tuple.
pub fn find_non_eulerian_binodes_with_differences<
    EdgeData: BidirectedData + Eq,
    Graph: StaticEdgeCentricBigraph<EdgeData = EdgeData>,
>(
    graph: &Graph,
) -> Vec<(Graph::NodeIndex, isize)> {
    let mut node_indices_and_differences = Vec::new();
    for node_index in graph.node_indices() {
        let difference = compute_eulerian_superfluous_out_biedges(graph, node_index);
        if difference != 0 {
            if graph.is_self_mirror_node(node_index) {
                node_indices_and_differences.push((node_index, 0));
            } else {
                node_indices_and_differences.push((node_index, difference));
            }
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

        let e1 = bigraph.add_edge(n1, n3, ());
        let _e2 = bigraph.add_edge(n4, n2, ());
        let e3 = bigraph.add_edge(n3, n5, ());
        let _e4 = bigraph.add_edge(n6, n4, ());
        let e5 = bigraph.add_edge(n5, n1, ());
        let _e6 = bigraph.add_edge(n2, n6, ());

        let euler_cycles = compute_minimum_bidirected_eulerian_cycle_decomposition(&bigraph);
        debug_assert_eq!(euler_cycles, vec![vec![e1, e3, e5]]);
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

        let e1 = bigraph.add_edge(n1, n3, ());
        let _e2 = bigraph.add_edge(n4, n2, ());
        let e3 = bigraph.add_edge(n3, n5, ());
        let _e4 = bigraph.add_edge(n6, n4, ());
        let e5 = bigraph.add_edge(n5, n1, ());
        let _e6 = bigraph.add_edge(n2, n6, ());
        let e7 = bigraph.add_edge(n1, n7, ());
        let _e8 = bigraph.add_edge(n8, n2, ());
        let e9 = bigraph.add_edge(n7, n9, ());
        let _e10 = bigraph.add_edge(n10, n8, ());
        let e11 = bigraph.add_edge(n9, n1, ());
        let _e12 = bigraph.add_edge(n2, n10, ());

        let euler_cycles = compute_minimum_bidirected_eulerian_cycle_decomposition(&bigraph);
        debug_assert_eq!(euler_cycles, vec![vec![e1, e3, e5, e7, e9, e11]]);
    }
}
