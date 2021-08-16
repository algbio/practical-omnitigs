use super::{MaximalMicrotigsAlgorithm, Microtigs};
use crate::macrotigs::macronodes::Macronodes;
use crate::restricted_reachability::{
    compute_inverse_restricted_forward_reachability, compute_restricted_backward_reachability,
    compute_restricted_forward_reachability,
};
use traitgraph::implementation::bit_vector_subgraph::BitVectorSubgraph;
use traitgraph::interface::subgraph::DecoratingSubgraph;
use traitgraph::interface::NavigableGraph;
use traitgraph::interface::StaticGraph;
use traitgraph::walks::VecEdgeWalk;
use traitsequence::interface::Sequence;

/// Compute the maximal microtigs of a strongly connected graph using hydrostructure-based queries.
pub struct StronglyConnectedHydrostructureBasedMaximalMicrotigs;

impl<Graph: StaticGraph> MaximalMicrotigsAlgorithm<Graph>
    for StronglyConnectedHydrostructureBasedMaximalMicrotigs
{
    fn compute_maximal_microtigs(
        graph: &Graph,
        macronodes: &Macronodes<Graph>,
    ) -> Microtigs<Graph> {
        let mut result = Vec::new();
        let mut macronodes_without_microtig_amount = 0;

        for macronode in macronodes.iter() {
            debug_assert!(!macronode.is_empty());
            let mut has_microtig = false;
            let center_in_node = *macronode.iter().next().unwrap();
            let center_out_node = *macronode.iter().last().unwrap();

            //////////////////////////////////////////
            //// Compute right-maximal microtigs. ////
            //////////////////////////////////////////

            // Choose any outgoing edge of the last node of the macronode center, and compute its
            // inverse restricted forwards reachability.
            let out_edge = graph
                .out_neighbors(center_out_node)
                .next()
                .expect(
                    "Found sink, but this algorithm requires the graph to be strongly connected",
                )
                .edge_id;
            let inverse_r_plus: BitVectorSubgraph<_> =
                compute_inverse_restricted_forward_reachability(graph, out_edge);

            // From the inverse restricted forward reachability, there are two cases:
            //  * It found exactly one incoming edge into the macronode center. Then that edge is a candidate for a separate central-micro omnitig.
            //  * It found more than one incoming edge into the macronode center. Then the only possible central-micro omnitig ends with out_edge.
            debug_assert!(inverse_r_plus.contains_node(center_in_node) || center_in_node == center_out_node, "Inverse restricted forward reachability did not find the first node of the macronode center, even though the macronode center has more than one node.");
            let in_edge_candidates: Vec<_> = inverse_r_plus
                .in_neighbors(center_in_node)
                .map(|n| n.edge_id)
                .collect();
            debug_assert!(!in_edge_candidates.is_empty(), "The inverse restricted reachability did not find any incoming edge into the macronode center.");

            // Compute maximal microtig based on out_edge.
            let mlmo1 = extend_left_micro_omnitig(graph, out_edge);
            let mrmo1 = if let Some(mlmo1) = &mlmo1 {
                extend_right_micro_omnitig(
                    graph,
                    *mlmo1
                        .get(mlmo1.len() - macronode.len() - 1)
                        .expect("Found left-micro omnitig, but it appears to be too short."),
                )
            } else {
                None
            };
            debug_assert_eq!(
                mlmo1.is_some(),
                mrmo1.is_some(),
                "Found either only a left-micro or only a right-micro omnitig."
            );
            if let (Some(mlmo1), Some(mrmo1)) = (&mlmo1, &mrmo1) {
                debug_assert_eq!(
                    mlmo1
                        .get(mlmo1.len() - macronode.len() - 1)
                        .expect("Left-micro omnitig misses macronode."),
                    mrmo1.iter().next().expect("Right-micro omnitig is empty."),
                    "Left-micro and right-micro omnitigs do not match."
                );
                debug_assert_eq!(
                    mlmo1.iter().last().expect("Left-micro omnitig is empty."),
                    mrmo1
                        .get(macronode.len())
                        .expect("Right-micro omnitig misses macronode."),
                    "Left-micro and right-micro omnitigs do not match."
                );
                result.push(
                    mlmo1
                        .iter()
                        .take(mlmo1.len() - macronode.len() - 1)
                        .chain(mrmo1.iter())
                        .copied()
                        .collect::<Vec<_>>(),
                );
                has_microtig = true;
            }

            if in_edge_candidates.len() == 1 {
                let in_edge = *in_edge_candidates.first().unwrap();
                let skip = if let Some(mrmo1) = &mrmo1 {
                    in_edge == *mrmo1.iter().next().unwrap()
                } else {
                    // TODO I am not sure if we actually cannot skip in this case. Not skipping is not wrong though, might just be slower.
                    false
                };

                if !skip {
                    let mrmo2 = extend_right_micro_omnitig(graph, in_edge);
                    let mlmo2 = if let Some(mrmo2) = &mrmo2 {
                        extend_left_micro_omnitig(
                            graph,
                            *mrmo2.get(macronode.len()).expect(
                                "Found right-micro omnitig, but it appears to be too short.",
                            ),
                        )
                    } else {
                        None
                    };
                    debug_assert_eq!(
                        mlmo2.is_some(),
                        mrmo2.is_some(),
                        "Found either only a left-micro or only a right-micro omnitig."
                    );
                    if let (Some(mlmo2), Some(mrmo2)) = (mlmo2, mrmo2) {
                        debug_assert_eq!(
                            mlmo2
                                .get(mlmo2.len() - macronode.len() - 1)
                                .expect("Left-micro omnitig misses macronode."),
                            mrmo2.get(0).expect("Right-micro omnitig is empty."),
                            "Left-micro and right-micro omnitigs do not match: {:?} -> {:?} via {:?}", mlmo2, mrmo2, macronode
                        );
                        debug_assert_eq!(
                            mlmo2.iter().last().expect("Left-micro omnitig is empty."),
                            mrmo2
                                .get(macronode.len())
                                .expect("Right-micro omnitig misses macronode."),
                            "Left-micro and right-micro omnitigs do not match."
                        );
                        result.push(
                            mlmo2
                                .iter()
                                .take(mlmo2.len() - macronode.len() - 1)
                                .chain(mrmo2.iter())
                                .copied()
                                .collect::<Vec<_>>(),
                        );
                        has_microtig = true;
                    }
                }
            }

            if !has_microtig {
                macronodes_without_microtig_amount += 1;
            }
        }

        Microtigs::from_vec_with_statistics(result, macronodes_without_microtig_amount)
    }
}

/// Returns the right-micro omnitig starting with the given edge.
/// In case the returned omnitig would be trivial, we return `None`, since then it does not have its macronode center as proper subwalk, so it is not a right-micro omnitig.
///
/// To handle the bivalent paths between maximal microtigs correctly, a right-micro omnitig is defined to end in the last edge of a bivalent path if it ends in a bivalent path.
fn extend_right_micro_omnitig<Graph: StaticGraph>(
    graph: &Graph,
    first_edge: Graph::EdgeIndex,
) -> Option<VecEdgeWalk<Graph>> {
    let r_minus: BitVectorSubgraph<_> = compute_restricted_backward_reachability(graph, first_edge);
    let mut rmo = vec![first_edge];
    let mut current_node = graph.edge_endpoints(first_edge).to_node;
    // Set to true if our omnitig is non-trivial. We return None otherwise.
    let mut non_trivial = false;

    // Move forward along the unique edge in r_minus, until it is not unique anymore or we encountered a join node.
    loop {
        let mut out_neighbors = r_minus.out_neighbors(current_node);
        let next_neighbor = out_neighbors.next().expect("No out neighbor in r_minus. This means that there either is a bug or the original graph was not strongly connected.");

        // Check if there is a unique next node.
        if out_neighbors.next().is_some() {
            break;
        }

        // Check if we are advancing past a split node.
        if !non_trivial && graph.out_degree(current_node) > 1 {
            non_trivial = true;
        }

        // Advance microtig.
        rmo.push(next_neighbor.edge_id);
        current_node = next_neighbor.node_id;

        // Check for bivalent edge.
        if graph.in_degree(current_node) > 1 {
            break;
        }
    }

    if non_trivial {
        Some(rmo)
    } else {
        None
    }
}

/// Returns the left-micro omnitig ending with the given edge.
/// In case the returned omnitig would be trivial, we return `None`, since then it does not have its macronode center as proper subwalk, so it is not a left-micro omnitig.
///
/// To handle the bivalent paths between maximal microtigs correctly, a left-micro omnitig is defined to end in the last edge of a bivalent path if it ends in a bivalent path.
fn extend_left_micro_omnitig<Graph: StaticGraph>(
    graph: &Graph,
    first_edge: Graph::EdgeIndex,
) -> Option<VecEdgeWalk<Graph>> {
    let r_plus: BitVectorSubgraph<_> = compute_restricted_forward_reachability(graph, first_edge);
    let mut lmo = vec![first_edge];
    // Since we want to end in the last edge of a bivalent path, we collect our edges in this vector first,
    // and only flush them whenever we found a join edge.
    let mut potential_bivalent_prefix = Vec::new();
    let mut current_node = graph.edge_endpoints(first_edge).from_node;
    // Set to true if our omnitig is non-trivial. We return None otherwise.
    let mut non_trivial = false;

    // Move forward along the unique edge in r_minus, until it is not unique anymore or we encountered a join node.
    loop {
        let mut in_neighbors = r_plus.in_neighbors(current_node);
        let next_neighbor = in_neighbors.next().expect("No out neighbor in r_plus. This means that there either is a bug or the original graph was not strongly connected.");

        // Check if there is a unique next node.
        if in_neighbors.next().is_some() {
            break;
        }

        // Check if we are advancing past a join node.
        let flush = graph.in_degree(current_node) > 1;
        if flush {
            non_trivial = true;
        }

        // Advance microtig.
        potential_bivalent_prefix.push(next_neighbor.edge_id);
        current_node = next_neighbor.node_id;

        // If we just added a join edge, our prefix is not bivalent.
        if flush {
            lmo.append(&mut potential_bivalent_prefix);
        }

        // Check for bivalent edge.
        if graph.out_degree(current_node) > 1 {
            break;
        }
    }

    if non_trivial {
        lmo.reverse();
        Some(lmo)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::MutableGraphContainer;
    use crate::macrotigs::macronodes::strongly_connected_macronode_algorithm::StronglyConnectedMacronodes;
    use crate::macrotigs::microtigs::strongly_connected_hydrostructure_based_maximal_microtig_algorithm::StronglyConnectedHydrostructureBasedMaximalMicrotigs;
    use crate::macrotigs::macronodes::MacronodeAlgorithm;
    use crate::macrotigs::microtigs::MaximalMicrotigsAlgorithm;
    use crate::macrotigs::microtigs::Microtigs;
    use traitgraph::interface::WalkableGraph;

    #[test]
    fn test_compute_maximal_microtigs_one_secluded() {
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
        debug_assert_eq!(
            maximal_microtigs,
            Microtigs::from(vec![graph.create_edge_walk(&[e6, e0, e1, e2, e10])])
        );
    }

    #[test]
    fn test_compute_maximal_microtigs_one_with_cross_bivalent() {
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
        debug_assert_eq!(
            maximal_microtigs,
            Microtigs::from(vec![
                graph.create_edge_walk(&[e6, e0, e1, e2, e10, e22, e23, e24]),
                graph.create_edge_walk(&[e24, e25]),
                graph.create_edge_walk(&[e29, e28, e27, e26, e6]),
            ])
        );
    }

    #[test]
    fn test_compute_maximal_microtigs_two() {
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
        debug_assert_eq!(
            maximal_microtigs,
            Microtigs::from(vec![
                graph.create_edge_walk(&[e9, e0, e2, e4]),
                graph.create_edge_walk(&[e10, e0, e1, e3]),
                graph.create_edge_walk(&[e5, e7, e9]),
                graph.create_edge_walk(&[e3, e5]),
                graph.create_edge_walk(&[e6, e8, e10]),
                graph.create_edge_walk(&[e4, e6]),
            ])
        );
    }
}
