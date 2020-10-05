use crate::hydrostructure::incremental_hydrostructure::NodeBridgeLikeIncrementalHydrostructure;
use crate::macrotigs::macrotigs::Macrotigs;
use crate::omnitigs::default_trivial_omnitigs::is_edge_in_maximal_trivial_omnitig_heart;
use bitvector::BitVector;
use traitgraph::algo::components::is_strong_bridge;
use traitgraph::index::GraphIndex;
use traitgraph::interface::StaticGraph;
use traitgraph::walks::{EdgeWalk, VecEdgeWalk};
use traitgraph::walks::{NodeWalk, VecNodeWalk};
use traitsequence::interface::Sequence;

fn check_safety<'graph, 'walk, Graph: StaticGraph>(
    graph: &'graph Graph,
    incremental_hydrostructure: &NodeBridgeLikeIncrementalHydrostructure<'graph, 'walk, Graph>,
) -> bool
where
    Graph::NodeIndex: 'static,
{
    let current_node_walk: VecNodeWalk<_> = incremental_hydrostructure
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
        for edge in macrotig.iter() {
            used_edges.insert(edge.as_usize());
        }
        safe_walks.extend(
            compute_maximal_node_covering_node_visible_one_circular_safe_subwalks(graph, &macrotig),
        );
    }

    for edge in graph.edge_indices() {
        if used_edges.contains(edge.as_usize())
            || !is_edge_in_maximal_trivial_omnitig_heart(graph, edge)
        {
            continue;
        }

        let trivial_omnitig: VecEdgeWalk<Graph> =
            EdgeWalk::compute_univocal_extension((vec![edge]).as_slice(), graph);
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

    safe_walks
}

/// Computes the maximal subwalks that are safe under the node-covering node-visible 1-circular walk model.
fn compute_maximal_node_covering_node_visible_one_circular_safe_subwalks<Graph: StaticGraph>(
    graph: &Graph,
    walk: &VecEdgeWalk<Graph>,
) -> Vec<VecNodeWalk<Graph>>
where
    Graph::NodeIndex: 'static,
{
    assert!(
        !walk.is_empty(),
        "Cannot compute safe subwalks of an empty walk."
    );
    let mut safe_walks = Vec::new();

    if walk.len() == 1 {
        let edge = walk.first().expect("Walk is empty.");
        if is_strong_bridge(graph, *edge) {
            safe_walks.push(NodeWalk::compute_univocal_extension(
                &[*edge]
                    .clone_as_node_walk::<VecNodeWalk<_>>(graph)
                    .expect("Walk cannot be represented as node walk."),
                graph,
            ));
        }

        return safe_walks;
    }

    // This reallocates memory every loop. It might make sense to allow to reuse the same structures for multiple walks.
    let mut incremental_hydrostructure =
        NodeBridgeLikeIncrementalHydrostructure::compute_and_set_fingers_left(graph, walk);

    // On a macrotig we can assume that length-2 walks are always bridge-like.
    while incremental_hydrostructure.can_increment_right_finger()
        || !check_safety(graph, &incremental_hydrostructure)
    {
        if check_safety(graph, &incremental_hydrostructure) {
            incremental_hydrostructure.increment_right_finger();

            if !check_safety(graph, &incremental_hydrostructure) {
                let safe_walk = incremental_hydrostructure.current_walk();
                let safe_walk = &safe_walk[0..safe_walk.len() - 1];
                safe_walks.push(NodeWalk::compute_univocal_extension(
                    &safe_walk
                        .clone_as_node_walk::<VecNodeWalk<_>>(graph)
                        .expect("Walk cannot be represented as node walk."),
                    graph,
                ));
            }
        } else if incremental_hydrostructure.can_increment_left_finger() {
            incremental_hydrostructure.increment_left_finger();
        } else {
            incremental_hydrostructure.increment_right_finger();
            incremental_hydrostructure.increment_left_finger();

            let edge = incremental_hydrostructure
                .current_walk()
                .first()
                .expect("Current walk is empty.");
            if !check_safety(graph, &incremental_hydrostructure) && is_strong_bridge(graph, *edge) {
                safe_walks.push(NodeWalk::compute_univocal_extension(
                    &[*edge]
                        .clone_as_node_walk::<VecNodeWalk<_>>(graph)
                        .expect("Walk cannot be represented as node walk."),
                    graph,
                ));
            }
        }
    }

    safe_walks.push(NodeWalk::compute_univocal_extension(
        &incremental_hydrostructure
            .current_walk()
            .clone_as_node_walk::<VecNodeWalk<_>>(graph)
            .expect("Walk cannot be represented as node walk."),
        graph,
    ));

    safe_walks
}
