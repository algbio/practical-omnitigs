use crate::hydrostructure::incremental_hydrostructure::NodeBridgeLikeIncrementalHydrostructure;
use crate::macrotigs::macrotigs::Macrotigs;
use crate::omnitigs::default_trivial_omnitigs::is_edge_in_maximal_trivial_omnitig_heart;
use crate::walks::{EdgeOmnitigLikeExt, NodeOmnitigLikeExt};
use bitvector::BitVector;
use traitgraph::index::GraphIndex;
use traitgraph::interface::StaticGraph;
use traitgraph::walks::{EdgeWalk, VecEdgeWalk};
use traitgraph::walks::{NodeWalk, VecNodeWalk};
use traitgraph_algo::components::is_strong_bridge;
use traitsequence::interface::Sequence;

fn check_safety<'graph, 'walk, Graph: StaticGraph>(
    graph: &'graph Graph,
    incremental_hydrostructure: &NodeBridgeLikeIncrementalHydrostructure<'graph, 'walk, Graph>,
) -> bool
where
    Graph::NodeIndex: 'static,
{
    let current_node_walk: VecNodeWalk<Graph> = incremental_hydrostructure
        .current_walk()
        .clone_as_node_walk(graph)
        .expect("Walk cannot be represented as node walk.");
    incremental_hydrostructure.is_safe()
        || current_node_walk
            .compute_trivial_heart_node_len(graph)
            .unwrap_or(0)
            > 0
    //<[Graph::NodeIndex] as NodeWalk<Graph>>::compute_trivial_heart_node_len(current_node_walk, graph).unwrap_or(0) > 0
}

/// Computes the maximal walks that are safe under the node-covering node-visible 1-circular walk model.
/// Does not return single nodes.
pub fn compute_maximal_node_covering_node_visible_one_circular_safe_walks<Graph: StaticGraph>(
    graph: &Graph,
    macrotigs: &Macrotigs<Graph>,
) -> Vec<VecNodeWalk<Graph>>
where
    Graph::NodeIndex: 'static,
{
    let mut safe_walks = Vec::new();
    let mut used_edges = BitVector::new(graph.edge_count());
    for macrotig in macrotigs.iter() {
        let univocal_extension =
            EdgeOmnitigLikeExt::compute_univocal_extension::<VecEdgeWalk<Graph>>(macrotig, graph);
        for edge in univocal_extension.iter() {
            used_edges.insert(edge.as_usize());
        }
        safe_walks.extend(
            compute_maximal_node_covering_node_visible_one_circular_safe_subwalks(graph, macrotig),
        );
    }

    for edge in graph.edge_indices() {
        if used_edges.contains(edge.as_usize())
            || !is_edge_in_maximal_trivial_omnitig_heart(graph, edge)
        {
            continue;
        }

        let trivial_omnitig: VecEdgeWalk<Graph> =
            EdgeOmnitigLikeExt::compute_univocal_extension((vec![edge]).as_slice(), graph);
        for edge in trivial_omnitig.iter() {
            used_edges.insert(edge.as_usize());
        }

        safe_walks.extend(
            compute_maximal_node_covering_node_visible_one_circular_safe_subwalks(
                graph,
                &trivial_omnitig,
            ),
        );
    }

    // This algorithm might produce non-maximal safe walks, therefore we need to filter them.
    // TODO this is asymptotically definitely not optimal
    let mut delete_indices = vec![false; safe_walks.len()];
    for (i, walk) in safe_walks.iter().enumerate() {
        let mut delete = false;
        for (j, superwalk) in safe_walks.iter().enumerate().filter(|(j, _)| i != *j) {
            if (walk == superwalk && i < j)
                || NodeWalk::<Graph, [Graph::NodeIndex]>::is_proper_subwalk_of(walk, superwalk)
            {
                delete = true;
                break;
            }
        }
        if delete {
            delete_indices[i] = true;
        }
    }

    let mut filtered_safe_walks = Vec::new();
    for (delete, walk) in delete_indices.iter().zip(safe_walks.into_iter()) {
        if !delete {
            filtered_safe_walks.push(walk);
        }
    }

    filtered_safe_walks
}

/// Computes the maximal subwalks that are safe under the node-covering node-visible 1-circular walk model.
fn compute_maximal_node_covering_node_visible_one_circular_safe_subwalks<Graph: StaticGraph>(
    graph: &Graph,
    walk: &[Graph::EdgeIndex],
) -> Vec<VecNodeWalk<Graph>>
where
    Graph::NodeIndex: 'static,
{
    debug_assert!(
        !walk.is_empty(),
        "Cannot compute safe subwalks of an empty walk."
    );
    let mut safe_walks = Vec::new();

    if walk.len() == 1 {
        let edge = walk.first().expect("Walk is empty.");
        if is_strong_bridge(graph, *edge) {
            safe_walks.push(NodeOmnitigLikeExt::compute_univocal_extension(
                &[*edge]
                    .clone_as_node_walk::<VecNodeWalk<Graph>>(graph)
                    .expect("Walk cannot be represented as node walk."),
                graph,
            ));
        }

        return safe_walks;
    }

    let mut incremental_hydrostructure =
        NodeBridgeLikeIncrementalHydrostructure::compute_and_set_fingers_left(graph, walk);

    let mut first_loop = true;
    while incremental_hydrostructure.can_increment_right_finger()
        || !check_safety(graph, &incremental_hydrostructure)
    {
        if check_safety(graph, &incremental_hydrostructure) {
            incremental_hydrostructure.increment_right_finger();

            if !check_safety(graph, &incremental_hydrostructure) {
                let safe_walk = incremental_hydrostructure.current_walk();
                let safe_walk = &safe_walk[0..safe_walk.len() - 1];
                safe_walks.push(NodeOmnitigLikeExt::compute_univocal_extension(
                    &safe_walk
                        .clone_as_node_walk::<VecNodeWalk<Graph>>(graph)
                        .expect("Walk cannot be represented as node walk."),
                    graph,
                ));
            }
        } else if incremental_hydrostructure.can_increment_left_finger() {
            incremental_hydrostructure.increment_left_finger();
        } else if incremental_hydrostructure.can_increment_right_finger() {
            // If the first two arcs are not safe together, then we need to check if the first arc is safe separately.
            if first_loop {
                let edge = incremental_hydrostructure
                    .current_walk()
                    .first()
                    .expect("Current walk is empty.");
                if !check_safety(graph, &incremental_hydrostructure)
                    && is_strong_bridge(graph, *edge)
                {
                    safe_walks.push(NodeOmnitigLikeExt::compute_univocal_extension(
                        &[*edge]
                            .clone_as_node_walk::<VecNodeWalk<Graph>>(graph)
                            .expect("Walk cannot be represented as node walk."),
                        graph,
                    ));
                }
            }

            incremental_hydrostructure.increment_right_finger();
            incremental_hydrostructure.increment_left_finger();

            // If two arcs are not safe together, and the following two-arc subwalk is also not safe,
            // then we need to check if the arc shared between the two is safe on its own.
            let edge = incremental_hydrostructure
                .current_walk()
                .first()
                .expect("Current walk is empty.");
            if !check_safety(graph, &incremental_hydrostructure) && is_strong_bridge(graph, *edge) {
                safe_walks.push(NodeOmnitigLikeExt::compute_univocal_extension(
                    &[*edge]
                        .clone_as_node_walk::<VecNodeWalk<Graph>>(graph)
                        .expect("Walk cannot be represented as node walk."),
                    graph,
                ));
            }
        } else {
            // If we are at the last two arcs and they are not safe, then check if the last arc is safe on its own.
            // Then, return the safe arcs directly, and do not add the last two arcs as safe as would be done otherwise.
            let edge = incremental_hydrostructure
                .current_walk()
                .last()
                .expect("Current walk is empty.");
            if is_strong_bridge(graph, *edge) {
                safe_walks.push(NodeOmnitigLikeExt::compute_univocal_extension(
                    &[*edge]
                        .clone_as_node_walk::<VecNodeWalk<Graph>>(graph)
                        .expect("Walk cannot be represented as node walk."),
                    graph,
                ));
            }

            return safe_walks;
        }

        first_loop = false;
    }

    // When the loop aborts we have a suffix of the walk that is maximal safe.
    // Therefore, we add it to the safe walks.
    safe_walks.push(NodeOmnitigLikeExt::compute_univocal_extension(
        &incremental_hydrostructure
            .current_walk()
            .clone_as_node_walk::<VecNodeWalk<Graph>>(graph)
            .expect("Walk cannot be represented as node walk."),
        graph,
    ));

    safe_walks
}

#[cfg(test)]
mod tests {
    use crate::macrotigs::macrotigs::Macrotigs;
    use crate::node_covering_node_visible_one_circular_safe::compute_maximal_node_covering_node_visible_one_circular_safe_walks;
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::MutableGraphContainer;
    use traitgraph::interface::WalkableGraph;

    #[test]
    fn test_compute_node_centric_omnitigs_simple() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let n4 = graph.add_node(());
        let n5 = graph.add_node(());
        let n6 = graph.add_node(());
        let n7 = graph.add_node(());
        let n8 = graph.add_node(());
        let n9 = graph.add_node(());
        let n10 = graph.add_node(());
        let n11 = graph.add_node(());
        let n12 = graph.add_node(());
        let n13 = graph.add_node(());
        let n14 = graph.add_node(());
        let n15 = graph.add_node(());
        let n16 = graph.add_node(());
        let n17 = graph.add_node(());
        let n18 = graph.add_node(());
        let n19 = graph.add_node(());
        let n20 = graph.add_node(());

        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n2, n3, ());
        let _e3 = graph.add_edge(n2, n4, ());
        let _e4 = graph.add_edge(n2, n5, ());
        let _e5 = graph.add_edge(n2, n6, ());
        let e6 = graph.add_edge(n7, n0, ()); // Comes from all except n11.
        let _e7 = graph.add_edge(n8, n0, ());
        let _e8 = graph.add_edge(n9, n0, ());
        let _e9 = graph.add_edge(n10, n0, ());
        let e10 = graph.add_edge(n3, n11, ()); // Goes to all except n7.
        let _e11 = graph.add_edge(n3, n12, ());
        let _e12 = graph.add_edge(n4, n13, ());
        let _e13 = graph.add_edge(n4, n14, ());
        let _e14 = graph.add_edge(n17, n8, ());
        let _e15 = graph.add_edge(n17, n9, ());
        let _e16 = graph.add_edge(n17, n10, ());
        let _e17 = graph.add_edge(n12, n18, ());
        let _e18 = graph.add_edge(n13, n18, ());
        let _e19 = graph.add_edge(n14, n18, ());
        let _e20 = graph.add_edge(n5, n18, ());
        let _e21 = graph.add_edge(n6, n18, ());
        let e22 = graph.add_edge(n11, n15, ());
        let e23 = graph.add_edge(n15, n16, ());
        let e24 = graph.add_edge(n16, n17, ());
        let e25 = graph.add_edge(n17, n17, ());
        let e26 = graph.add_edge(n20, n7, ());
        let e27 = graph.add_edge(n19, n20, ());
        let e28 = graph.add_edge(n18, n19, ());
        let e29 = graph.add_edge(n18, n18, ());

        let maximal_macrotigs = Macrotigs::compute(&graph);
        debug_assert_eq!(
            maximal_macrotigs,
            Macrotigs::from(vec![graph.create_edge_walk(&[
                e29, e28, e27, e26, e6, e0, e1, e2, e10, e22, e23, e24, e25
            ]),])
        );

        let maximal_node_centric_omnitigs =
            compute_maximal_node_covering_node_visible_one_circular_safe_walks(
                &graph,
                &maximal_macrotigs,
            );
        debug_assert_eq!(
            maximal_node_centric_omnitigs,
            vec![
                Vec::from([n18, n19, n20, n7, n0, n1, n2, n3, n11, n15, n16, n17]),
                graph.create_node_walk(&[n0, n1, n2, n5, n18]),
                graph.create_node_walk(&[n0, n1, n2, n6, n18]),
                graph.create_node_walk(&[n17, n8, n0, n1, n2]),
                graph.create_node_walk(&[n17, n9, n0, n1, n2]),
                graph.create_node_walk(&[n17, n10, n0, n1, n2]),
                graph.create_node_walk(&[n0, n1, n2, n3, n12, n18]),
                graph.create_node_walk(&[n0, n1, n2, n4, n13, n18]),
                graph.create_node_walk(&[n0, n1, n2, n4, n14, n18]),
            ]
        );
    }

    #[test]
    fn test_minimal_trivial_walk_with_unsafe_heart() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());

        let _e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n2, n3, ());
        let e3 = graph.add_edge(n3, n1, ());
        let e4 = graph.add_edge(n3, n0, ());
        let e5 = graph.add_edge(n0, n2, ());

        let maximal_macrotigs = Macrotigs::compute(&graph);
        debug_assert_eq!(
            maximal_macrotigs,
            Macrotigs::from(vec![graph.create_edge_walk(&[e3, e1, e2, e4, e5])])
        );

        let maximal_node_centric_omnitigs =
            compute_maximal_node_covering_node_visible_one_circular_safe_walks(
                &graph,
                &maximal_macrotigs,
            );
        debug_assert_eq!(
            maximal_node_centric_omnitigs,
            vec![Vec::from([n1, n2, n3, n0]),]
        );
    }

    #[test]
    fn test_trivial_node_heart() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());

        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n0, ());
        let e2 = graph.add_edge(n0, n0, ());

        let maximal_macrotigs = Macrotigs::compute(&graph);
        debug_assert_eq!(
            maximal_macrotigs,
            Macrotigs::from(vec![
                graph.create_edge_walk(&[e1, e2]),
                graph.create_edge_walk(&[e2, e0, e1]),
            ])
        );

        let maximal_node_centric_omnitigs =
            compute_maximal_node_covering_node_visible_one_circular_safe_walks(
                &graph,
                &maximal_macrotigs,
            );
        debug_assert_eq!(maximal_node_centric_omnitigs, vec![Vec::from([n0, n1, n0])]);
    }

    #[test]
    fn test_trivial_sea_cloud_node() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());

        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n2, n3, ());
        let e3 = graph.add_edge(n3, n0, ());
        let e4 = graph.add_edge(n1, n0, ());
        let e5 = graph.add_edge(n3, n2, ());

        let maximal_macrotigs = Macrotigs::compute(&graph);
        debug_assert_eq!(
            maximal_macrotigs,
            Macrotigs::from(vec![
                graph.create_edge_walk(&[e5, e2, e3, e0, e4]),
                graph.create_edge_walk(&[e4, e0, e1, e2, e5]),
            ])
        );

        let maximal_node_centric_omnitigs =
            compute_maximal_node_covering_node_visible_one_circular_safe_walks(
                &graph,
                &maximal_macrotigs,
            );
        debug_assert_eq!(
            maximal_node_centric_omnitigs,
            vec![
                Vec::from([n2, n3, n0, n1]),
                graph.create_node_walk(&[n0, n1, n2, n3]),
            ],
        );
    }
}
