use super::{Macrotigs, MaximalMacrotigsAlgorithm};
use crate::macrotigs::microtigs::Microtigs;
use traitgraph::index::GraphIndex;
use traitgraph::index::OptionalGraphIndex;
use traitgraph::interface::StaticGraph;
use traitgraph_algo::traversal::univocal_traversal::is_edge_self_bivalent;
use traitsequence::interface::Sequence;

/// Compute the maximal microtigs of a strongly connected graph using hydrostructure-based queries.
pub struct DefaultMacrotigLinkAlgorithm;

impl<Graph: StaticGraph> MaximalMacrotigsAlgorithm<Graph> for DefaultMacrotigLinkAlgorithm {
    fn compute_maximal_macrotigs(graph: &Graph, microtigs: &Microtigs<Graph>) -> Macrotigs<Graph> {
        let mut result = Vec::new();
        let mut outgoing_microtigs = vec![Graph::OptionalEdgeIndex::new_none(); graph.edge_count()];
        let mut incoming_microtigs = vec![Graph::OptionalEdgeIndex::new_none(); graph.edge_count()];
        let mut used_microtigs = vec![false; microtigs.len()];

        // Prepare mapping from edges to incoming and outgoing microtigs.
        for (microtig_index, microtig) in microtigs.iter().enumerate() {
            debug_assert!(!microtig.is_empty(), "Found empty microtig.");
            let first_edge = microtig.first().unwrap();
            let last_edge = microtig.last().unwrap();

            debug_assert!(outgoing_microtigs[first_edge.as_usize()].is_none());
            debug_assert!(incoming_microtigs[last_edge.as_usize()].is_none());
            outgoing_microtigs[first_edge.as_usize()] = Some(microtig_index).into();
            incoming_microtigs[last_edge.as_usize()] = Some(microtig_index).into();
        }

        // Combine microtigs.
        for (microtig_index, microtig) in microtigs.iter().enumerate() {
            if !used_microtigs[microtig_index] {
                //println!("Using microtig {}", microtig_index);

                used_microtigs[microtig_index] = true;
                let mut macrotig = vec![microtig_index];
                let mut current_microtig = microtig;

                // Extend to the left.
                loop {
                    let first_edge = current_microtig.first().unwrap();

                    // Do not extend over self bivalent edges as required by the definition of macrotigs.
                    if is_edge_self_bivalent(graph, *first_edge) {
                        break;
                    }

                    if let Some(incoming_microtig_index) =
                        incoming_microtigs[first_edge.as_usize()].as_usize()
                    {
                        /*println!(
                            "Found incoming microtig {} at edge {}",
                            incoming_microtig_index,
                            first_edge.as_usize()
                        );*/

                        // Append incoming microtig and make it new current.
                        macrotig.push(incoming_microtig_index);
                        current_microtig = &microtigs[incoming_microtig_index];
                        debug_assert!(
                            !used_microtigs[incoming_microtig_index],
                            "Trying to use a microtig for two maximal macrotigs."
                        );
                        used_microtigs[incoming_microtig_index] = true;
                    } else {
                        break;
                    }
                }

                macrotig.reverse();
                current_microtig = microtig;

                // Extend to the right.
                loop {
                    let last_edge = current_microtig.last().unwrap();

                    // Do not extend over self bivalent edges as required by the definition of macrotigs.
                    if is_edge_self_bivalent(graph, *last_edge) {
                        break;
                    }

                    if let Some(outgoing_microtig_index) =
                        outgoing_microtigs[last_edge.as_usize()].as_usize()
                    {
                        /*println!(
                            "Found outgoing microtig {} at edge {}",
                            outgoing_microtig_index,
                            last_edge.as_usize()
                        );*/

                        // Append incoming microtig and make it new current.
                        macrotig.push(outgoing_microtig_index);
                        current_microtig = &microtigs[outgoing_microtig_index];
                        debug_assert!(
                            !used_microtigs[outgoing_microtig_index],
                            "Trying to use a microtig for two maximal macrotigs."
                        );
                        used_microtigs[outgoing_microtig_index] = true;
                    } else {
                        break;
                    }
                }

                // Combine microtigs into macrotigs.
                let mut macrotig_walk = vec![microtigs[macrotig[0]][0]];
                for microtig_index in macrotig {
                    let microtig = &microtigs[microtig_index];
                    macrotig_walk.extend(microtig.iter().skip(1));
                }
                result.push(macrotig_walk);
            }
        }

        Macrotigs::from(result)
    }
}

#[cfg(test)]
mod tests {
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::{MutableGraphContainer, WalkableGraph};
    use crate::macrotigs::macronodes::strongly_connected_macronode_algorithm::StronglyConnectedMacronodes;
    use crate::macrotigs::microtigs::strongly_connected_hydrostructure_based_maximal_microtig_algorithm::StronglyConnectedHydrostructureBasedMaximalMicrotigs;
    use crate::macrotigs::macronodes::MacronodeAlgorithm;
    use crate::macrotigs::microtigs::MaximalMicrotigsAlgorithm;
    use crate::macrotigs::macrotigs::default_macrotig_link_algorithm::DefaultMacrotigLinkAlgorithm;
    use crate::macrotigs::macrotigs::{MaximalMacrotigsAlgorithm, Macrotigs};

    #[test]
    fn test_compute_maximal_macrotigs_one_secluded() {
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
        let _e14 = graph.add_edge(n11, n8, ());
        let _e15 = graph.add_edge(n11, n9, ());
        let _e16 = graph.add_edge(n11, n10, ());
        let _e17 = graph.add_edge(n12, n7, ());
        let _e18 = graph.add_edge(n13, n7, ());
        let _e19 = graph.add_edge(n14, n7, ());
        let _e20 = graph.add_edge(n5, n7, ());
        let _e21 = graph.add_edge(n6, n7, ());

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
            Macrotigs::from(vec![graph.create_edge_walk(&[e6, e0, e1, e2, e10])])
        );
    }

    #[test]
    fn test_compute_maximal_macrotigs_one_with_cross_bivalent() {
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
    }

    #[test]
    fn test_compute_maximal_macrotigs_two() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let n4 = graph.add_node(());
        let n5 = graph.add_node(());
        let n6 = graph.add_node(());
        let n7 = graph.add_node(());

        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n1, n3, ());
        let e3 = graph.add_edge(n2, n6, ());
        let e4 = graph.add_edge(n3, n7, ());
        let e5 = graph.add_edge(n6, n6, ());
        let e6 = graph.add_edge(n7, n7, ());
        let e7 = graph.add_edge(n6, n4, ());
        let e8 = graph.add_edge(n7, n5, ());
        let e9 = graph.add_edge(n4, n0, ());
        let e10 = graph.add_edge(n5, n0, ());

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
            Macrotigs::from(vec![
                graph.create_edge_walk(&[e5, e7, e9, e0, e2, e4, e6]),
                graph.create_edge_walk(&[e6, e8, e10, e0, e1, e3, e5]),
            ])
        );
    }

    #[test]
    fn test_cycle_with_central_bypass() {
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
    }
}
