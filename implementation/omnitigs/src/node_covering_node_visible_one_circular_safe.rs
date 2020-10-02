use crate::hydrostructure::incremental_hydrostructure::NodeBridgeLikeIncrementalHydrostructure;
use crate::macrotigs::macrotigs::Macrotigs;
use traitgraph::interface::StaticGraph;
use traitgraph::walks::EdgeWalk;
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

/// Computes the maximal non-trivial walks that are safe under the node-covering node-visible 1-circular walk model.
pub fn compute_maximal_non_trivial_node_covering_node_visible_one_circular_safe_walks<
    Graph: StaticGraph,
>(
    graph: &Graph,
    macrotigs: &Macrotigs<Graph>,
) -> Vec<VecNodeWalk<Graph>>
where
    Graph::NodeIndex: 'static,
{
    let mut safe_walks = Vec::new();

    for macrotig in macrotigs.iter() {
        assert!(
            macrotig.len() >= 2,
            "Macrotigs have a length of at least two edges."
        );

        // This reallocates memory every loop. It might make sense to allow to reuse the same structures for multiple walks.
        let mut incremental_hydrostructure =
            NodeBridgeLikeIncrementalHydrostructure::compute_and_set_fingers_left(graph, macrotig);

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
                todo!("Handle len 1 walks");
            }
        }

        safe_walks.push(NodeWalk::compute_univocal_extension(
            &incremental_hydrostructure
                .current_walk()
                .clone_as_node_walk::<VecNodeWalk<_>>(graph)
                .expect("Walk cannot be represented as node walk."),
            graph,
        ));
    }

    safe_walks
}
