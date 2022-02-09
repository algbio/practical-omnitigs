use crate::omnitigs::univocal_extension_algorithms::{
    NonSccUnivocalExtensionStrategy, SccUnivocalExtensionStrategy,
};
use crate::omnitigs::{Omnitig, Omnitigs, TrivialOmnitigAlgorithm, UnivocalExtensionAlgorithm};
use bitvector::BitVector;
use std::marker::PhantomData;
use traitgraph::index::GraphIndex;
use traitgraph::interface::{NodeOrEdge, StaticGraph};
use traitgraph::walks::VecEdgeWalk;
use traitgraph_algo::traversal::univocal_traversal::UnivocalIterator;
use traitsequence::interface::Sequence;

/// An algorithm to extract trivial omnitigs.
pub struct DefaultTrivialOmnitigAlgorithm<UnivocalExtensionStrategy> {
    _univocal_extension_strategy: PhantomData<UnivocalExtensionStrategy>,
}

/// An algorithm to extract trivial omnitigs form a strongly connected graph.
pub type SccTrivialOmnitigAlgorithm = DefaultTrivialOmnitigAlgorithm<SccUnivocalExtensionStrategy>;

/// An algorithm to extract trivial omnitigs form a not strongly connected graph.
/// This runs slightly slower than the counterpart for strongly connected graphs, especially for long univocal extensions.
pub type NonSccTrivialOmnitigAlgorithm =
    DefaultTrivialOmnitigAlgorithm<NonSccUnivocalExtensionStrategy>;

/// Returns true if the edge is in a trivial omnitig heart.
pub fn is_edge_in_maximal_trivial_omnitig_heart<Graph: StaticGraph>(
    graph: &Graph,
    edge: Graph::EdgeIndex,
) -> bool {
    let edge_endpoints = graph.edge_endpoints(edge);

    let unitig_start_node =
        UnivocalIterator::new_backward_without_start(graph, NodeOrEdge::Edge(edge))
            .filter_map(|n_or_e| match n_or_e {
                NodeOrEdge::Node(node) => {
                    if !graph.is_biunivocal_node(node) || node == edge_endpoints.to_node {
                        Some(node)
                    } else {
                        None
                    }
                }
                _ => None,
            })
            .next()
            .unwrap();
    // True if the edge is part of a cycle that can be traversed R-univocally indefinitely.
    let reverse_univocal_cycle =
        unitig_start_node == edge_endpoints.to_node && graph.in_degree(unitig_start_node) == 1;

    let unitig_end_node =
        UnivocalIterator::new_forward_without_start(graph, NodeOrEdge::Edge(edge))
            .filter_map(|n_or_e| match n_or_e {
                NodeOrEdge::Node(node) => {
                    if !graph.is_biunivocal_node(node) || node == edge_endpoints.from_node {
                        Some(node)
                    } else {
                        None
                    }
                }
                _ => None,
            })
            .next()
            .unwrap();
    // True if the edge is part of a cycle that can be traversed univocally indefinitely.
    let forward_univocal_cycle =
        unitig_end_node == edge_endpoints.from_node && graph.out_degree(unitig_end_node) == 1;
    // True if the edge is part of a biunivocal cycle. That is a cycle that is disconnected to all other possibly existing parts of the graph.
    let univocal_cycle = forward_univocal_cycle && reverse_univocal_cycle;

    // True if the edge is part of a trivial heart in a non-circular weakly connected component of the graph.
    // This means that its surrounding unitig cannot be entered via a univocal walk, and not entered in reverse via an R-univocal walk.
    let heart_in_non_circular_graph = (graph.out_degree(unitig_start_node) >= 2
        || graph.in_degree(unitig_start_node) == 0)
        && (graph.in_degree(unitig_end_node) >= 2 || graph.out_degree(unitig_end_node) == 0);

    univocal_cycle || heart_in_non_circular_graph
}

impl<
        Graph: StaticGraph,
        UnivocalExtensionStrategy: UnivocalExtensionAlgorithm<Graph, VecEdgeWalk<Graph>>,
    > TrivialOmnitigAlgorithm<Graph> for DefaultTrivialOmnitigAlgorithm<UnivocalExtensionStrategy>
{
    type UnivocalExtensionStrategy = UnivocalExtensionStrategy;

    fn compute_maximal_trivial_omnitigs(
        graph: &Graph,
        mut omnitigs: Omnitigs<Graph>,
    ) -> Omnitigs<Graph> {
        debug!("Marking used edges");
        let mut used_edges = BitVector::new(graph.edge_count());
        for omnitig in omnitigs.iter() {
            for edge in omnitig.iter() {
                used_edges.insert(edge.as_usize());
            }
        }

        debug!("Extend {} unused edges", graph.edge_count());
        for edge in graph.edge_indices() {
            if used_edges.contains(edge.as_usize())
                || !is_edge_in_maximal_trivial_omnitig_heart(graph, edge)
            {
                continue;
            }

            let trivial_omnitig: VecEdgeWalk<Graph> =
                Self::UnivocalExtensionStrategy::compute_univocal_extension(graph, &[edge]);
            for edge in trivial_omnitig.iter() {
                used_edges.insert(edge.as_usize());
            }
            let last_split_edge = trivial_omnitig
                .iter()
                .enumerate()
                .filter(|(_, e)| graph.is_split_edge(**e))
                .map(|(i, _)| i)
                .last()
                .unwrap_or(0);
            let first_join_edge = trivial_omnitig
                .iter()
                .enumerate()
                .filter(|(_, e)| graph.is_join_edge(**e))
                .map(|(i, _)| i)
                .next()
                .unwrap_or(trivial_omnitig.len() - 1);

            omnitigs.push(Omnitig::new(
                trivial_omnitig,
                last_split_edge,
                first_join_edge,
            ));
        }

        debug_assert_eq!(used_edges.len(), graph.edge_count());

        omnitigs
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
    use crate::omnitigs::{MacrotigBasedNonTrivialOmnitigAlgorithm, Omnitigs, TrivialOmnitigAlgorithm, Omnitig};
    use crate::omnitigs::default_trivial_omnitigs::SccTrivialOmnitigAlgorithm;

    #[test]
    fn test_compute_omnitigs_simple() {
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
        let e3 = graph.add_edge(n2, n4, ());
        let e4 = graph.add_edge(n2, n5, ());
        let e5 = graph.add_edge(n2, n6, ());
        let e6 = graph.add_edge(n7, n0, ()); // Comes from all except n11.
        let e7 = graph.add_edge(n8, n0, ());
        let e8 = graph.add_edge(n9, n0, ());
        let e9 = graph.add_edge(n10, n0, ());
        let e10 = graph.add_edge(n3, n11, ()); // Goes to all except n7.
        let e11 = graph.add_edge(n3, n12, ());
        let e12 = graph.add_edge(n4, n13, ());
        let e13 = graph.add_edge(n4, n14, ());
        let e14 = graph.add_edge(n17, n8, ());
        let e15 = graph.add_edge(n17, n9, ());
        let e16 = graph.add_edge(n17, n10, ());
        let e17 = graph.add_edge(n12, n18, ());
        let e18 = graph.add_edge(n13, n18, ());
        let e19 = graph.add_edge(n14, n18, ());
        let e20 = graph.add_edge(n5, n18, ());
        let e21 = graph.add_edge(n6, n18, ());
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

        let maximal_omnitigs = SccTrivialOmnitigAlgorithm::compute_maximal_trivial_omnitigs(
            &graph,
            maximal_non_trivial_omnitigs,
        );
        debug_assert_eq!(
            maximal_omnitigs,
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
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e4, e20]), 2, 3),
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e5, e21]), 2, 3),
                Omnitig::new(graph.create_edge_walk(&[e14, e7, e0, e1]), 0, 1),
                Omnitig::new(graph.create_edge_walk(&[e15, e8, e0, e1]), 0, 1),
                Omnitig::new(graph.create_edge_walk(&[e16, e9, e0, e1]), 0, 1),
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e2, e11, e17]), 3, 4),
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e3, e12, e18]), 3, 4),
                Omnitig::new(graph.create_edge_walk(&[e0, e1, e3, e13, e19]), 3, 4),
            ])
        );
    }

    #[test]
    fn test_compute_trivial_omnitigs_cycle() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n2, n0, ());

        let trivial_omnitigs = Omnitigs::compute_trivial_only(&graph);
        debug_assert_eq!(
            trivial_omnitigs,
            Omnitigs::from(vec![Omnitig::new(
                graph.create_edge_walk(&[e1, e2, e0, e1, e2]),
                0,
                4
            )])
        );
    }

    #[test]
    fn test_compute_trivial_omnitigs_trap_cycle() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let e0 = graph.add_edge(n3, n0, ());
        let e1 = graph.add_edge(n0, n1, ());
        let e2 = graph.add_edge(n1, n2, ());
        let e3 = graph.add_edge(n2, n0, ());

        let trivial_omnitigs = Omnitigs::compute_trivial_only_non_scc(&graph);
        debug_assert_eq!(
            trivial_omnitigs,
            Omnitigs::from(vec![Omnitig::new(
                graph.create_edge_walk(&[e0, e1, e2, e3]),
                0,
                0
            )])
        );
    }

    #[test]
    fn test_compute_trivial_omnitigs_reverse_trap_cycle() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let e0 = graph.add_edge(n0, n3, ());
        let e1 = graph.add_edge(n0, n1, ());
        let e2 = graph.add_edge(n1, n2, ());
        let e3 = graph.add_edge(n2, n0, ());

        let trivial_omnitigs = Omnitigs::compute_trivial_only_non_scc(&graph);
        debug_assert_eq!(
            trivial_omnitigs,
            Omnitigs::from(vec![Omnitig::new(
                graph.create_edge_walk(&[e1, e2, e3, e0]),
                3,
                3
            )])
        );
    }

    #[test]
    fn test_compute_trivial_omnitigs_trap_cycle_inverse_id_order() {
        let mut graph = petgraph_impl::new();
        let n3 = graph.add_node(());
        let n2 = graph.add_node(());
        let n1 = graph.add_node(());
        let n0 = graph.add_node(());
        let e3 = graph.add_edge(n2, n0, ());
        let e2 = graph.add_edge(n1, n2, ());
        let e1 = graph.add_edge(n0, n1, ());
        let e0 = graph.add_edge(n3, n0, ());

        let trivial_omnitigs = Omnitigs::compute_trivial_only_non_scc(&graph);
        debug_assert_eq!(
            trivial_omnitigs,
            Omnitigs::from(vec![Omnitig::new(
                graph.create_edge_walk(&[e0, e1, e2, e3]),
                0,
                0
            )])
        );
    }

    #[test]
    fn test_compute_trivial_omnitigs_reverse_trap_cycle_inverse_id_order() {
        let mut graph = petgraph_impl::new();
        let n3 = graph.add_node(());
        let n2 = graph.add_node(());
        let n1 = graph.add_node(());
        let n0 = graph.add_node(());
        let e3 = graph.add_edge(n2, n0, ());
        let e2 = graph.add_edge(n1, n2, ());
        let e1 = graph.add_edge(n0, n1, ());
        let e0 = graph.add_edge(n0, n3, ());

        let trivial_omnitigs = Omnitigs::compute_trivial_only_non_scc(&graph);
        debug_assert_eq!(
            trivial_omnitigs,
            Omnitigs::from(vec![Omnitig::new(
                graph.create_edge_walk(&[e1, e2, e3, e0]),
                3,
                3
            )])
        );
    }

    #[test]
    fn test_compute_trivial_omnitigs_path() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());

        let trivial_omnitigs = Omnitigs::compute_trivial_only_non_scc(&graph);
        debug_assert_eq!(
            trivial_omnitigs,
            Omnitigs::from(vec![Omnitig::new(graph.create_edge_walk(&[e0, e1]), 0, 1)])
        );
    }
}
