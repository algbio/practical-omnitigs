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
use traitgraph::walks::EdgeWalk;
use traitgraph::walks::NodeWalk;
use traitgraph::walks::VecEdgeWalk;

/// Compute the maximal microtigs of a strongly connected graph using hydrostructure-based queries.
pub struct StronglyConnectedHydrostructureBasedMaximalMicrotigs;

impl<Graph: StaticGraph> MaximalMicrotigsAlgorithm<Graph>
    for StronglyConnectedHydrostructureBasedMaximalMicrotigs
/*where
Graph::EdgeIndex: 'static,*/
{
    fn compute_maximal_microtigs(
        graph: &Graph,
        macronodes: &Macronodes<Graph>,
    ) -> Microtigs<Graph> {
        let mut result = Vec::new();

        for macronode in macronodes {
            assert!(!macronode.is_empty());
            let center_in_node = macronode.iter().next().unwrap();
            let center_out_node = macronode.iter().last().unwrap();

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
            assert!(inverse_r_plus.contains_node(center_in_node) || center_in_node == center_out_node, "Inverse restricted forward reachability did not find the first node of the macronode center, even though the macronode center has more than one node.");
            let in_edge_candidates: Vec<_> = inverse_r_plus
                .in_neighbors(center_in_node)
                .map(|n| n.edge_id)
                .collect();
            assert!(!in_edge_candidates.is_empty(), "The inverse restricted reachability did not find any incoming edge into the macronode center.");

            // Compute maximal microtig based on out_edge.
            let mlmo1 = extend_left_micro_omnitig(graph, out_edge);
            let mrmo1 = if let Some(mlmo1) = &mlmo1 {
                extend_right_micro_omnitig(
                    graph,
                    mlmo1
                        .iter()
                        .nth(mlmo1.len() - macronode.len() - 1)
                        .expect("Found left-micro omnitig, but it appears to be too short."),
                )
            } else {
                None
            };
            assert_eq!(
                mlmo1.is_some(),
                mrmo1.is_some(),
                "Found either only a left-micro or only a right-micro omnitig."
            );
            if let (Some(mlmo1), Some(mrmo1)) = (&mlmo1, &mrmo1) {
                assert_eq!(
                    mlmo1
                        .iter()
                        .nth(mlmo1.len() - macronode.len() - 1)
                        .expect("Left-micro omnitig misses macronode."),
                    mrmo1.iter().next().expect("Right-micro omnitig is empty."),
                    "Left-micro and right-micro omnitigs do not match."
                );
                assert_eq!(
                    mlmo1.iter().last().expect("Left-micro omnitig is empty."),
                    mrmo1
                        .iter()
                        .nth(macronode.len())
                        .expect("Right-micro omnitig misses macronode."),
                    "Left-micro and right-micro omnitigs do not match."
                );
                result.push(VecEdgeWalk::<Graph>::from(
                    mlmo1
                        .iter()
                        .take(mlmo1.len() - macronode.len() - 1)
                        .chain(mrmo1.iter())
                        .collect::<Vec<_>>(),
                ));
            }

            if in_edge_candidates.len() == 1 {
                let in_edge = *in_edge_candidates.first().unwrap();
                let skip = if let Some(mrmo1) = &mrmo1 {
                    in_edge == mrmo1.iter().next().unwrap()
                } else {
                    // TODO I am not sure if we actually cannot skip in this case. Not skipping is not wrong though, might just be slower.
                    false
                };

                if !skip {
                    let mrmo2 = extend_right_micro_omnitig(graph, in_edge);
                    let mlmo2 = if let Some(mrmo2) = &mrmo2 {
                        extend_left_micro_omnitig(
                            graph,
                            mrmo2.iter().nth(macronode.len()).expect(
                                "Found right-micro omnitig, but it appears to be too short.",
                            ),
                        )
                    } else {
                        None
                    };
                    assert_eq!(
                        mlmo2.is_some(),
                        mrmo2.is_some(),
                        "Found either only a left-micro or only a right-micro omnitig."
                    );
                    if let (Some(mlmo2), Some(mrmo2)) = (mlmo2, mrmo2) {
                        assert_eq!(
                            mlmo2
                                .iter()
                                .nth(mlmo2.len() - macronode.len() - 1)
                                .expect("Left-micro omnitig misses macronode."),
                            mrmo2.iter().next().expect("Right-micro omnitig is empty."),
                            "Left-micro and right-micro omnitigs do not match."
                        );
                        assert_eq!(
                            mlmo2.iter().last().expect("Left-micro omnitig is empty."),
                            mrmo2
                                .iter()
                                .nth(macronode.len())
                                .expect("Right-micro omnitig misses macronode."),
                            "Left-micro and right-micro omnitigs do not match."
                        );
                        result.push(VecEdgeWalk::<Graph>::from(
                            mlmo2
                                .iter()
                                .take(mlmo2.len() - macronode.len() - 1)
                                .chain(mrmo2.iter())
                                .collect::<Vec<_>>(),
                        ));
                    }
                }
            }
        }

        Microtigs::new(result)
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
        Some(rmo.into())
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
        Some(lmo.into())
    } else {
        None
    }
}
