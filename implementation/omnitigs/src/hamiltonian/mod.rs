use crate::macrotigs::macrotigs::Macrotigs;
use crate::node_covering_node_visible_one_circular_safe::compute_maximal_node_covering_node_visible_one_circular_safe_walks;
use bitvector::BitVector;
use std::collections::{BTreeMap, BTreeSet};
use traitgraph::algo::components::{is_cycle, is_strongly_connected};
use traitgraph::index::GraphIndex;
use traitgraph::index::OptionalGraphIndex;
use traitgraph::interface::{DynamicGraph, ImmutableGraphContainer, StaticGraph};
use traitgraph::walks::VecNodeWalk;
use traitsequence::interface::Sequence;

enum WalkOverlap {
    Forward,
    Backward,
    Both,
    None,
    Invalid,
}

fn check_walk_overlap<Graph: StaticGraph>(
    walk: &VecNodeWalk<Graph>,
    other_walk: &VecNodeWalk<Graph>,
) -> WalkOverlap {
    // Check if the walks share any nodes.
    let nodes = walk.iter().copied().collect::<BTreeSet<_>>();
    let other_nodes = other_walk.iter().copied().collect::<BTreeSet<_>>();
    let intersection = nodes
        .intersection(&other_nodes)
        .copied()
        .collect::<BTreeSet<_>>();

    if intersection.is_empty() {
        return WalkOverlap::None;
    }

    // Check if the intersection can be a valid forward overlap (i.e. an overlap which is a prefix of one walk and a suffix of another).
    let forward_overlap = if intersection.contains(walk.last().unwrap())
        && intersection.contains(other_walk.first().unwrap())
    {
        let overlap_start = walk
            .iter()
            .enumerate()
            .filter(|(_, node)| *node == other_walk.first().unwrap())
            .map(|(i, _)| i)
            .next()
            .unwrap();
        let overlap_end = other_walk
            .iter()
            .enumerate()
            .filter(|(_, node)| *node == walk.last().unwrap())
            .map(|(i, _)| i)
            .next()
            .unwrap();
        let other_overlap = &other_walk[..=overlap_end];
        let overlap = &walk[overlap_start..];

        assert!(
            overlap_start != 0 && overlap_end != walk.len() - 1,
            "One walk is a subwalk of the other."
        );

        if overlap != other_overlap {
            // If the other walk ends inside the walk and the walk starts inside the other walk, but this does not create a valid overlap, then the overlap is invalid.
            // Observe that they then cannot overlap only in the first/last node.
            return WalkOverlap::Invalid;
        }

        Some((overlap_start, overlap_end))
    } else {
        None
    };

    // Check if the intersection can be a valid backward overlap (i.e. an overlap which is a prefix of one walk and a suffix of another).
    let backward_overlap = if intersection.contains(walk.first().unwrap())
        && intersection.contains(other_walk.last().unwrap())
    {
        let overlap_start = other_walk
            .iter()
            .enumerate()
            .filter(|(_, node)| *node == walk.first().unwrap())
            .map(|(i, _)| i)
            .next()
            .unwrap();
        let overlap_end = walk
            .iter()
            .enumerate()
            .filter(|(_, node)| *node == other_walk.last().unwrap())
            .map(|(i, _)| i)
            .next()
            .unwrap();
        let other_overlap = &other_walk[overlap_start..];
        let overlap = &walk[..=overlap_end];

        assert!(
            overlap_start != 0 && overlap_end != walk.len() - 1,
            "One walk is a subwalk of the other."
        );

        if overlap != other_overlap {
            // If the other walk ends inside the walk and the walk starts inside the other walk, but this does not create a valid overlap, then the overlap is invalid.
            // Observe that they then cannot overlap only in the first/last node.
            return WalkOverlap::Invalid;
        }

        Some((overlap_start, overlap_end))
    } else {
        None
    };

    if let (Some((_, forward_overlap_end)), Some((backward_overlap_start, _))) =
        (forward_overlap, backward_overlap)
    {
        let walk = VecNodeWalk::<Graph>::new(
            walk.iter()
                .chain(
                    other_walk
                        .iter()
                        .take(backward_overlap_start + 1)
                        .skip(forward_overlap_end + 1),
                )
                .copied()
                .collect::<Vec<_>>(),
        );
        if check_walk_invalid_self_overlap(&walk) {
            WalkOverlap::Invalid
        } else {
            WalkOverlap::Both
        }
    } else if let Some((_, forward_overlap_end)) = forward_overlap {
        let walk = VecNodeWalk::<Graph>::new(
            walk.iter()
                .chain(other_walk.iter().skip(forward_overlap_end + 1))
                .copied()
                .collect::<Vec<_>>(),
        );
        if check_walk_invalid_self_overlap(&walk) {
            WalkOverlap::Invalid
        } else {
            WalkOverlap::Forward
        }
    } else if let Some((_, backward_overlap_end)) = backward_overlap {
        let walk = VecNodeWalk::<Graph>::new(
            other_walk
                .iter()
                .chain(walk.iter().skip(backward_overlap_end + 1))
                .copied()
                .collect::<Vec<_>>(),
        );
        if check_walk_invalid_self_overlap(&walk) {
            WalkOverlap::Invalid
        } else {
            WalkOverlap::Backward
        }
    } else {
        // If there is an overlap, but it cannot be made into a valid overlap, then it is invalid.
        WalkOverlap::Invalid
    }
}

/// Returns true if a walk overlaps itself in an invalid way.
fn check_walk_invalid_self_overlap<Graph: StaticGraph>(walk: &VecNodeWalk<Graph>) -> bool {
    let mut intersection = BTreeSet::new();
    for (i, node) in walk.iter().copied().enumerate() {
        if walk[..i].contains(&node) || walk[i + 1..].contains(&node) {
            intersection.insert(node);
        }
    }

    if intersection.is_empty() {
        return false;
    }

    if !intersection.contains(walk.first().unwrap()) {
        return false;
    }

    intersection.remove(walk.first().unwrap());
    !intersection.is_empty()
}

fn do_walks_cover_all_nodes<Graph: ImmutableGraphContainer>(
    graph: &Graph,
    walks: Vec<VecNodeWalk<Graph>>,
) -> bool {
    let mut covered_nodes = BitVector::new(graph.node_count());
    for walk in &walks {
        for node in walk.iter() {
            covered_nodes.insert(node.as_usize());
        }
    }
    covered_nodes.len() == graph.node_count()
}

/// Returns a reduced graph that is hamiltonian if and only if the input graph is hamiltonian,
/// or returns `None` if the graph can be proven to not be hamiltonian (which equals to the reduction being impossible).
pub fn preprocess_hamiltonian_circuit<Graph: DynamicGraph + Default>(graph: &Graph) -> Option<Graph>
where
    Graph::NodeIndex: 'static,
    Graph::NodeData: Default,
    Graph::EdgeData: Default,
{
    if !is_strongly_connected(graph) {
        return None;
    }

    if is_cycle(graph) {
        let mut result = Graph::default();
        result.add_node(Default::default());
        return Some(result);
    }

    if graph.node_count() <= 1 {
        let mut result = Graph::default();
        result.add_node(Default::default());
        return Some(result);
    }

    let safe_walks = compute_maximal_node_covering_node_visible_one_circular_safe_walks(
        graph,
        &Macrotigs::compute(graph),
    );

    if safe_walks
        .iter()
        .any(|walk| check_walk_invalid_self_overlap(walk))
    {
        return None;
    }

    // Map nodes to walks that contain them.
    let mut node_walk_map = vec![Vec::new(); graph.node_count()];
    for (i, walk) in safe_walks.iter().enumerate() {
        for node in walk.iter() {
            node_walk_map[node.as_usize()].push(i);
        }
    }
    let node_walk_map = node_walk_map;

    // Check if there are invalid overlaps or mergeable walks.
    let mut forward_mergeable_walks = BTreeMap::new();
    let mut backward_mergeable_walks = BTreeMap::new();
    let mut harmless_pairs = BTreeMap::new();
    for node in graph.node_indices() {
        // Check all walks that contain this node.
        for (skip_index, &walk_index) in node_walk_map[node.as_usize()].iter().enumerate() {
            for &other_walk_index in node_walk_map[node.as_usize()].iter().skip(skip_index + 1) {
                // If a pair of walks was already identified as mergeable or harmless, skip it.
                if forward_mergeable_walks.get(&walk_index) == Some(&other_walk_index)
                    || backward_mergeable_walks.get(&walk_index) == Some(&other_walk_index)
                    || harmless_pairs.get(&walk_index) == Some(&other_walk_index)
                {
                    continue;
                }

                let walk = &safe_walks[walk_index];
                let other_walk = &safe_walks[other_walk_index];

                match check_walk_overlap(walk, other_walk) {
                    WalkOverlap::Forward => {
                        let old_value =
                            forward_mergeable_walks.insert(walk_index, other_walk_index);
                        assert!(old_value.is_none(), "Found two walks that can be merged after a walk, which is a contradiction to the mergeability of either.");
                    }
                    WalkOverlap::Backward => {
                        let old_value =
                            backward_mergeable_walks.insert(other_walk_index, walk_index);
                        assert!(old_value.is_none(), "Found two walks that can be merged before a walk, which is a contradiction to the mergeability of either.");
                    }
                    WalkOverlap::Both => {
                        let old_value =
                            forward_mergeable_walks.insert(walk_index, other_walk_index);
                        assert!(old_value.is_none(), "Found two walks that can be merged after a walk, which is a contradiction to the mergeability of either.");
                        let old_value =
                            backward_mergeable_walks.insert(other_walk_index, walk_index);
                        assert!(old_value.is_none(), "Found two walks that can be merged before a walk, which is a contradiction to the mergeability of either.");
                    }
                    WalkOverlap::None => {
                        harmless_pairs.insert(walk_index, other_walk_index);
                        harmless_pairs.insert(other_walk_index, walk_index);
                    }
                    WalkOverlap::Invalid => {
                        return None;
                    }
                };
            }
        }
    }

    // Here we know that all pairs of walks are either harmless or mergeable.
    // First, we merge all mergeable walks.
    let mut merged_walks = Vec::new();
    for (i, walk) in safe_walks.iter().enumerate() {
        // Find beginning of merge.
        let mut first_walk_index = i;
        let mut circular = false;
        while let Some(predecessor_index) = backward_mergeable_walks.get(&first_walk_index).copied()
        {
            first_walk_index = predecessor_index;
            if first_walk_index == i {
                // Found cycle
                circular = true;
                break;
            }
        }

        // Collect walk indices and remove mappings.
        let mut walk_indices = vec![first_walk_index];
        while let Some(successor_index) = forward_mergeable_walks
            .get(&walk_indices.last().unwrap())
            .copied()
        {
            forward_mergeable_walks.remove(&walk_indices.last().unwrap());
            backward_mergeable_walks.remove(&successor_index);
            walk_indices.push(successor_index);
        }

        if circular {
            assert!(
                forward_mergeable_walks.is_empty()
                    && backward_mergeable_walks.is_empty()
                    && walk_indices.first() == walk_indices.last()
                    && walk_indices.len() == safe_walks.len() + 1,
                "Walks can be merged circularly, but there are other pairs left"
            );
            // If there is a (merged) circular safe walk, and this covers all nodes, then the graph is hamiltonian.
            // If the walk does not cover all nodes, then the graph is not hamiltonian.
            if do_walks_cover_all_nodes(graph, safe_walks) {
                let mut result = Graph::default();
                result.add_node(Default::default());
                return Some(result);
            } else {
                return None;
            }
        } else if walk_indices.len() == 1 {
            merged_walks.push(walk.clone())
        } else {
            let mut merged_walk = walk.clone();
            for walk_index in walk_indices.iter().copied().skip(1) {
                merged_walk = VecNodeWalk::<Graph>::new(
                    merged_walk
                        .forward_merge_iter_assume_mergeable(&safe_walks[walk_index])
                        .copied()
                        .collect::<Vec<_>>(),
                );
            }
            merged_walks.push(merged_walk);
        }
    }

    // Remove all inner nodes and their incident arcs of merged walks and instead insert and arc from the first to the last node of the merged walk.
    // Do this by copying the graph.
    let mut result = Graph::default();
    let mut node_id_map = vec![Graph::OptionalNodeIndex::new_none(); graph.node_count()];
    let mut processed_nodes = BitVector::new(graph.node_count());

    // Add merged walks as edges into the result graph.
    for walk in &merged_walks {
        for node in walk.iter() {
            assert!(
                !processed_nodes.contains(node.as_usize()),
                "A node is part of two separate not mergeable walks."
            );
            processed_nodes.insert(node.as_usize());
        }

        let result_first_node = result.add_node(Default::default());
        let result_last_node = result.add_node(Default::default());
        result.add_edge(result_first_node, result_last_node, Default::default());
        node_id_map[walk.first().unwrap().as_usize()] = Some(result_first_node).into();
        node_id_map[walk.last().unwrap().as_usize()] = Some(result_last_node).into();
    }

    // Add nodes that are not part of any merged walk to the result graph.
    for node in graph.node_indices() {
        if !processed_nodes.contains(node.as_usize()) {
            let result_node = result.add_node(Default::default());
            node_id_map[node.as_usize()] = Some(result_node).into();
        }
    }

    // Add edges between non-deleted nodes.
    for node in graph.node_indices() {
        if let Some(result_node) = node_id_map[node.as_usize()].into() {
            for neighbor in graph.out_neighbors(node) {
                if let Some(result_neighbor_node) = node_id_map[neighbor.node_id.as_usize()].into()
                {
                    result.add_edge(result_node, result_neighbor_node, Default::default());
                }
            }
        }
    }

    Some(result)
}

#[cfg(test)]
mod tests {
    use crate::hamiltonian::preprocess_hamiltonian_circuit;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use traitgraph::algo::predefined_graphs::{
        create_random_graph, create_random_hamiltonian_graph,
    };
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::MutableGraphContainer;

    #[test]
    fn test_random_hamiltonian_graphs() {
        let mut graph = petgraph_impl::new::<(), ()>();
        let mut random = StdRng::seed_from_u64(0);

        for _ in 0..10 {
            graph.clear();
            create_random_hamiltonian_graph(&mut graph, 20, 1.0, &mut random);
            let reduced = preprocess_hamiltonian_circuit(&graph);
            assert!(reduced.is_some());
        }
    }

    #[test]
    fn test_random_graphs() {
        let mut graph = petgraph_impl::new::<(), ()>();
        let mut random = StdRng::seed_from_u64(0);

        let mut result = Vec::new();
        for _ in 0..10 {
            graph.clear();
            create_random_graph(&mut graph, 20, 1.0, &mut random);
            let reduced = preprocess_hamiltonian_circuit(&graph);
            result.push(reduced.is_some());
        }

        // This assertion is dependent on the random generator used and might fail if its implementation changes.
        assert_eq!(
            result,
            vec![true, true, true, true, true, false, false, true, true, true]
        );
    }
}