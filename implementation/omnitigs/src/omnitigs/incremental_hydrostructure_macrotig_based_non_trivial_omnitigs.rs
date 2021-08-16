use crate::hydrostructure::incremental_hydrostructure::BridgeLikeIncrementalHydrostructure;
use crate::macrotigs::macrotigs::Macrotigs;
use crate::omnitigs::{MacrotigBasedNonTrivialOmnitigAlgorithm, Omnitig, Omnitigs};
use traitgraph::interface::StaticGraph;
use traitsequence::interface::Sequence;

/// A macrotig-based non-trivial omnitig algorithm that uses the incremental hydrostructure.
pub struct IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm;

impl<Graph: StaticGraph> MacrotigBasedNonTrivialOmnitigAlgorithm<Graph>
    for IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm
{
    fn compute_maximal_non_trivial_omnitigs(
        graph: &Graph,
        macrotigs: &Macrotigs<Graph>,
    ) -> Omnitigs<Graph> {
        let mut omnitigs = Vec::new();
        let mut omnitigs_per_macrotig = Vec::new();

        for macrotig in macrotigs.iter() {
            debug_assert!(
                macrotig.len() >= 2,
                "Macrotigs have a length of at least two edges."
            );

            // This reallocates memory every loop. It might make sense to allow to reuse the same structures for multiple walks.
            let mut incremental_hydrostructure =
                BridgeLikeIncrementalHydrostructure::compute_and_set_fingers_left(graph, macrotig);
            let mut omnitigs_per_macrotig_current = 0;

            // On a macrotig we can assume that length-2 walks are always bridge-like.
            while incremental_hydrostructure.can_increment_right_finger()
                || !incremental_hydrostructure.is_safe()
            {
                if incremental_hydrostructure.is_safe() {
                    incremental_hydrostructure.increment_right_finger();

                    if !incremental_hydrostructure.is_safe() {
                        let omnitig = incremental_hydrostructure.current_walk();
                        let omnitig = &omnitig[0..omnitig.len() - 1];
                        omnitigs.push(Omnitig::compute_from_non_trivial_heart_superwalk(
                            graph, omnitig,
                        ));
                        omnitigs_per_macrotig_current += 1;
                    }
                } else {
                    incremental_hydrostructure.increment_left_finger();
                }
            }

            omnitigs.push(Omnitig::compute_from_non_trivial_heart_superwalk(
                graph,
                incremental_hydrostructure.current_walk(),
            ));
            omnitigs_per_macrotig_current += 1;

            omnitigs_per_macrotig.push(omnitigs_per_macrotig_current);
        }

        Omnitigs::new(omnitigs, omnitigs_per_macrotig)
    }
}

#[cfg(test)]
mod tests {
    use traitgraph::implementation::petgraph_impl;
    use crate::macrotigs::macronodes::strongly_connected_macronode_algorithm::StronglyConnectedMacronodes;
    use crate::macrotigs::microtigs::strongly_connected_hydrostructure_based_maximal_microtig_algorithm::StronglyConnectedHydrostructureBasedMaximalMicrotigs;
    use crate::macrotigs::macronodes::MacronodeAlgorithm;
    use crate::macrotigs::microtigs::MaximalMicrotigsAlgorithm;
    use crate::macrotigs::macrotigs::default_macrotig_link_algorithm::DefaultMacrotigLinkAlgorithm;
    use crate::macrotigs::macrotigs::{MaximalMacrotigsAlgorithm, Macrotigs};
    use traitgraph::interface::WalkableGraph;
    use traitgraph::interface::MutableGraphContainer;
    use crate::omnitigs::incremental_hydrostructure_macrotig_based_non_trivial_omnitigs::IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm;
    use crate::omnitigs::{MacrotigBasedNonTrivialOmnitigAlgorithm, Omnitigs, Omnitig};

    #[test]
    fn test_compute_non_trivial_omnitigs_simple() {
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

        let macronodes = StronglyConnectedMacronodes::compute_macronodes(&graph);
        let maximal_microtigs =
            StronglyConnectedHydrostructureBasedMaximalMicrotigs::compute_maximal_microtigs(
                &graph,
                &macronodes,
            );
        let maximal_macrotigs =
            DefaultMacrotigLinkAlgorithm::compute_maximal_macrotigs(&graph, &maximal_microtigs);
        debug_assert_eq!(
            maximal_macrotigs,
            Macrotigs::from(vec![graph.create_edge_walk(&[
                e29, e28, e27, e26, e6, e0, e1, e2, e10, e22, e23, e24, e25
            ]),])
        );

        let maximal_non_trivial_omnitigs = IncrementalHydrostructureMacrotigBasedNonTrivialOmnitigAlgorithm::compute_maximal_non_trivial_omnitigs(&graph, &maximal_macrotigs);
        debug_assert_eq!(
            maximal_non_trivial_omnitigs,
            Omnitigs::from(vec![
                Omnitig::new(
                    graph.create_edge_walk(&[e29, e28, e27, e26, e6, e0, e1]),
                    0,
                    1
                ),
                Omnitig::new(
                    graph.create_edge_walk(&[e28, e27, e26, e6, e0, e1, e2, e10, e22, e23, e24]),
                    3,
                    7
                ),
                Omnitig::new(
                    graph.create_edge_walk(&[e0, e1, e2, e10, e22, e23, e24, e25]),
                    6,
                    7
                ),
            ])
        );
    }
}
