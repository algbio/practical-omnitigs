use crate::hydrostructure::Hydrostructure;
use crate::restricted_reachability::{
    compute_incremental_restricted_backward_edge_reachability,
    compute_incremental_restricted_forward_edge_reachability,
};
use traitgraph::implementation::incremental_subgraph::IncrementalSubgraph;
use traitgraph::interface::subgraph::DecoratingSubgraph;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::EdgeWalk;
use traitgraph::walks::VecEdgeWalk;

/// The hydrostructure for a walk `W`.
/// This hydrostructure implementation is incremental, meaning that it is valid for any subwalk of `W`.
/// The subwalk can be adjusted using the left and right finger.
pub struct IncrementalHydrostructure<'a, Graph: GraphBase> {
    /// An incremental version of the set `R⁺(W)` for each split of the underlying walk.
    r_plus: IncrementalSubgraph<'a, Graph>,
    /// An incremental version of the set `R⁻(W)` for each join of the underlying walk.
    r_minus: IncrementalSubgraph<'a, Graph>,
    /// The walk the incremental hydrostructure corresponds to.
    walk: VecEdgeWalk<Graph>,
    /// The first edge of the walk for which the hydrostructure should be valid.
    left_finger: usize,
    /// The last edge of the walk for which the hydrostructure should be valid.
    right_finger: usize,
    /// The rightmost split edge in `[left_finger + 1 ... right_finger]`.
    rightmost_split: Option<usize>,
    /// The rightmost join edge in `[left_finger ... right_finger - 1]`.
    rightmost_join: Option<usize>,
}

impl<'a, Graph: 'a + StaticGraph> IncrementalHydrostructure<'a, Graph> {
    /// Initialise the hydrostructure of a walk `W` with given sets `R⁺(W)`, `R⁻(W)`.
    /// Panics if the given walk has less than two edges, since the hydrostructure is defined only for walks of at least two edges.
    pub fn new(
        r_plus: IncrementalSubgraph<'a, Graph>,
        r_minus: IncrementalSubgraph<'a, Graph>,
        walk: VecEdgeWalk<Graph>,
    ) -> Self {
        assert!(
            walk.len() >= 2,
            "The hydrostructure is defined only for walks of at least two edges."
        );
        let mut result = Self {
            left_finger: 0,
            right_finger: walk.len() - 1,
            rightmost_split: None,
            rightmost_join: None,
            r_plus,
            r_minus,
            walk,
        };
        result.update_rightmost_split_join();
        result
    }

    /// Compute the incremental hydrostructure of a walk.
    /// Panics if the given walk has less than two edges, since the hydrostructure is defined only for walks of at least two edges.
    pub fn compute(graph: &'a Graph, walk: VecEdgeWalk<Graph>) -> Self {
        assert!(
            walk.len() >= 2,
            "The hydrostructure is defined only for walks of at least two edges."
        );
        let r_plus = compute_incremental_restricted_forward_edge_reachability(graph, &walk);
        let r_minus = compute_incremental_restricted_backward_edge_reachability(graph, &walk);
        let mut result = Self {
            left_finger: 0,
            right_finger: walk.len() - 1,
            rightmost_split: None,
            rightmost_join: None,
            r_plus,
            r_minus,
            walk,
        };
        result.update_rightmost_split_join();
        result
    }

    fn update_rightmost_split_join(&mut self) {
        self.rightmost_split = self
            .walk
            .iter()
            .enumerate()
            .take(self.right_finger + 1)
            .skip(self.left_finger + 1)
            .rev()
            .filter(|(_, e)| self.r_plus.parent_graph().is_split_edge(*e))
            .map(|(n, _)| n)
            .next();
        self.rightmost_join = self
            .walk
            .iter()
            .enumerate()
            .take(self.right_finger)
            .skip(self.left_finger)
            .rev()
            .filter(|(_, e)| self.r_plus.parent_graph().is_join_edge(*e))
            .map(|(n, _)| n)
            .next();
    }

    /// Set the left finger to the specified value.
    /// Panics if the given index is not a valid index in the underlying walk of this incremental hydrostructure, or if it is not left of the right finger.
    pub fn set_left_finger(&mut self, left_finger: usize) {
        assert!(left_finger < self.walk.len() && left_finger < self.right_finger);
        self.left_finger = left_finger;
        self.r_minus
            .set_current_step(self.walk.len() - 1 - left_finger);
        self.update_rightmost_split_join();
    }

    /// Returns the value of the left finger.
    pub fn left_finger(&self) -> usize {
        self.left_finger
    }

    /// Increments the left finger.
    /// Panics if the increment moves the left finger past the end of the walk.
    pub fn increment_left_finger(&mut self) {
        self.left_finger += 1;
        assert!(self.left_finger < self.walk.len() && self.left_finger < self.right_finger);
        self.r_minus
            .set_current_step(self.walk.len() - 1 - self.left_finger);

        if let Some(rightmost_split) = self.rightmost_split {
            if rightmost_split <= self.left_finger {
                self.rightmost_split = None;
            }
        }
        if let Some(rightmost_join) = self.rightmost_join {
            if rightmost_join < self.left_finger {
                self.rightmost_join = None;
            }
        }
    }

    /// Set the right finger to the specified value.
    /// Panics if the given index is not a valid index in the underlying walk of this incremental hydrostructure, or if it is not right of the left finger.
    pub fn set_right_finger(&mut self, right_finger: usize) {
        assert!(right_finger < self.walk.len() && self.left_finger < right_finger);
        self.right_finger = right_finger;
        self.r_plus.set_current_step(right_finger);
        self.update_rightmost_split_join();
    }

    /// Returns the value of the right finger.
    pub fn right_finger(&self) -> usize {
        self.right_finger
    }

    /// Increments the right finger.
    /// Panics if the increment moves the right finger past the end of the walk.
    pub fn increment_right_finger(&mut self) {
        if self
            .r_plus
            .parent_graph()
            .is_join_edge(self.walk[self.right_finger])
        {
            self.rightmost_join = Some(self.right_finger);
        }

        self.right_finger += 1;
        assert!(self.right_finger < self.walk.len());
        self.r_plus.set_current_step(self.right_finger);

        if self
            .r_plus
            .parent_graph()
            .is_split_edge(self.walk[self.right_finger])
        {
            self.rightmost_split = Some(self.right_finger);
        }
    }
}

impl<'a, Graph: 'a + StaticGraph> Hydrostructure<Graph::NodeIndex, Graph::EdgeIndex>
    for IncrementalHydrostructure<'a, Graph>
{
    fn is_node_r_plus(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        if self.r_plus.contains_node(node) && self.rightmost_split.is_some() {
            true
        } else {
            self.walk
                .iter()
                .take(self.right_finger)
                .skip(self.left_finger)
                .map(|e| self.r_plus.parent_graph().edge_endpoints(e).to_node)
                .any(|n| n == node)
        }
    }

    fn is_node_r_minus(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        if self.r_minus.contains_node(node) && self.rightmost_join.is_some() {
            true
        } else {
            self.walk
                .iter()
                .take(self.right_finger)
                .skip(self.left_finger)
                .map(|e| self.r_plus.parent_graph().edge_endpoints(e).to_node)
                .any(|n| n == node)
        }
    }

    fn is_edge_r_plus(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        if self.r_plus.contains_edge(edge) && self.rightmost_split.is_some() {
            true
        } else {
            self.walk
                .iter()
                .take(self.right_finger)
                .skip(self.left_finger)
                .any(|e| e == edge)
        }
    }

    fn is_edge_r_minus(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        if self.r_minus.contains_edge(edge) && self.rightmost_join.is_some() {
            true
        } else {
            self.walk
                .iter()
                .take(self.right_finger + 1)
                .skip(self.left_finger + 1)
                .any(|e| e == edge)
        }
    }

    fn is_bridge_like(&self) -> bool {
        unimplemented!(
            "This cannot be implemented efficiently for this structure without tracking the vapor"
        )
    }

    fn is_avertible(&self) -> bool {
        unimplemented!(
            "This cannot be implemented efficiently for this structure without tracking the vapor"
        )
    }
}

#[cfg(test)]
mod tests {
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::{MutableGraphContainer, WalkableGraph, ImmutableGraphContainer};
    use crate::macrotigs::macronodes::strongly_connected_macronode_algorithm::StronglyConnectedMacronodes;
    use crate::macrotigs::microtigs::strongly_connected_hydrostructure_based_maximal_microtig_algorithm::StronglyConnectedHydrostructureBasedMaximalMicrotigs;
    use crate::macrotigs::macronodes::MacronodeAlgorithm;
    use crate::macrotigs::microtigs::MaximalMicrotigsAlgorithm;
    use crate::macrotigs::macrotigs::default_macrotig_link_algorithm::DefaultMacrotigLinkAlgorithm;
    use crate::macrotigs::macrotigs::{MaximalMacrotigsAlgorithm, Macrotigs};
    use crate::hydrostructure::incremental_hydrostructure::IncrementalHydrostructure;
    use traitgraph::walks::{EdgeWalk, VecEdgeWalk};
    use crate::hydrostructure::static_hydrostructure::StaticHydrostructure;
    use crate::hydrostructure::Hydrostructure;

    #[test]
    fn test_incremental_hydrostructure_macrotig() {
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
        assert_eq!(
            maximal_macrotigs,
            Macrotigs::new(vec![graph.create_edge_walk(&[
                e29, e28, e27, e26, e6, e0, e1, e2, e10, e22, e23, e24, e25
            ]),])
        );

        for macrotig in maximal_macrotigs.iter() {
            let mut incremental_hydrostructure =
                IncrementalHydrostructure::compute(&graph, macrotig.clone());
            for subwalk_len in 2..macrotig.len() {
                //println!("Setting initial left finger for subwalk len {}", subwalk_len);
                incremental_hydrostructure.set_left_finger(0);
                //println!("Setting initial right finger for subwalk len {}", subwalk_len);
                incremental_hydrostructure.set_right_finger(subwalk_len - 1);
                for offset in 0..(macrotig.len() - subwalk_len) {
                    let subwalk: VecEdgeWalk<_> =
                        graph.create_edge_walk(&macrotig[offset..offset + subwalk_len]);
                    //println!("{:?}", subwalk);
                    //println!("has rightmost split/join: {:?}/{:?}; left/right finger: {}/{}", incremental_hydrostructure.rightmost_split, incremental_hydrostructure.rightmost_join, incremental_hydrostructure.left_finger, incremental_hydrostructure.right_finger);
                    let static_hydrostructure =
                        StaticHydrostructure::compute_with_bitvector_subgraph(
                            &graph,
                            subwalk.clone(),
                        );

                    if static_hydrostructure.is_bridge_like() {
                        for node in graph.node_indices() {
                            assert_eq!(
                                incremental_hydrostructure.is_node_r_plus(node),
                                static_hydrostructure.is_node_r_plus(node),
                                "Node {:?} r_plus incremental/static: {}/{}",
                                node,
                                incremental_hydrostructure.is_node_r_plus(node),
                                static_hydrostructure.is_node_r_plus(node)
                            );
                            assert_eq!(
                                incremental_hydrostructure.is_node_r_minus(node),
                                static_hydrostructure.is_node_r_minus(node),
                                "Node {:?} r_minus incremental/static: {}/{}",
                                node,
                                incremental_hydrostructure.is_node_r_minus(node),
                                static_hydrostructure.is_node_r_minus(node)
                            );
                        }
                        for edge in graph.edge_indices() {
                            assert_eq!(
                                incremental_hydrostructure.is_edge_r_plus(edge),
                                static_hydrostructure.is_edge_r_plus(edge),
                                "Edge {:?} r_plus incremental/static: {}/{}",
                                edge,
                                incremental_hydrostructure.is_edge_r_plus(edge),
                                static_hydrostructure.is_edge_r_plus(edge)
                            );
                            assert_eq!(
                                incremental_hydrostructure.is_edge_r_minus(edge),
                                static_hydrostructure.is_edge_r_minus(edge),
                                "Edge {:?} r_minus incremental/static: {}/{}",
                                edge,
                                incremental_hydrostructure.is_edge_r_minus(edge),
                                static_hydrostructure.is_edge_r_minus(edge)
                            );
                        }
                    }

                    if offset < (macrotig.len() - subwalk_len) - 1 {
                        incremental_hydrostructure.increment_right_finger();
                        incremental_hydrostructure.increment_left_finger();
                    }
                }
            }
        }
    }
}
