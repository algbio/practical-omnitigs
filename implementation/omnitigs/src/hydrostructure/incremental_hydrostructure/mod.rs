use crate::hydrostructure::incremental_hydrostructure::conjunctive_safety_tracker::ConjunctiveSafetyTracker;
use crate::hydrostructure::incremental_hydrostructure::node_centric_component_tracker::NodeCentricComponentTracker;
use crate::hydrostructure::Hydrostructure;
use crate::restricted_reachability::{
    compute_incremental_restricted_backward_edge_reachability,
    compute_incremental_restricted_forward_edge_reachability,
};
use traitgraph::implementation::incremental_subgraph::IncrementalSubgraph;
use traitgraph::interface::subgraph::DecoratingSubgraph;
use traitgraph::interface::{GraphBase, StaticGraph};
use vapor_is_path_tracker::VaporIsPathTracker;

/// A type that combines two safety trackers under conjunction.
pub mod conjunctive_safety_tracker;
/// A type that keeps counts of the nodes in the different hydrostructure components to dynamically determine if they contain nodes.
pub mod node_centric_component_tracker;
/// A type that keeps counts about nodes and edges in a subgraph to dynamically determine if the subgraph is a path.
pub mod vapor_is_path_tracker;

/// An incremental hydrostructure that checks if a subwalk is bridge-like.
pub type BridgeLikeIncrementalHydrostructure<'graph, 'walk, Graph> =
    IncrementalHydrostructure<'graph, 'walk, Graph, VaporIsPathTracker<'graph, Graph>>;

/// An incremental hydrostructure that checks if the node sequence (including the tail of the last arc and the head of the first arc)
/// of a subwalk is safe in the node-visible node-covering 1-closed walk model.
pub type NodeBridgeLikeIncrementalHydrostructure<'graph, 'walk, Graph> = IncrementalHydrostructure<
    'graph,
    'walk,
    Graph,
    ConjunctiveSafetyTracker<VaporIsPathTracker<'graph, Graph>, NodeCentricComponentTracker>,
>;

/// The hydrostructure for a walk `W`.
/// This hydrostructure implementation is incremental, meaning that it is valid for any subwalk of `W`.
/// The subwalk can be adjusted using the left and right finger.
pub struct IncrementalHydrostructure<
    'graph,
    'walk,
    Graph: GraphBase,
    SafetyTracker: IncrementalSafetyTracker<'graph, Graph>,
> {
    /// An incremental version of the set `R⁺(W)` for each split of the underlying walk.
    r_plus: IncrementalSubgraph<'graph, Graph>,
    /// An incremental version of the set `R⁻(W)` for each join of the underlying walk.
    r_minus: IncrementalSubgraph<'graph, Graph>,
    /// The walk the incremental hydrostructure corresponds to.
    walk: &'walk [Graph::EdgeIndex],
    /// The first edge of the walk for which the hydrostructure should be valid.
    left_finger: usize,
    /// The last edge of the walk for which the hydrostructure should be valid.
    right_finger: usize,
    /// The rightmost split edge in `[left_finger + 1 ... right_finger]`.
    rightmost_split: Option<usize>,
    /// The rightmost join edge in `[left_finger ... right_finger - 1]`.
    rightmost_join: Option<usize>,
    /// Track if the current subwalk is safe.
    safety_tracker: SafetyTracker,
}

impl<
        'graph,
        'walk,
        Graph: 'graph + StaticGraph,
        SafetyTracker: IncrementalSafetyTracker<'graph, Graph>,
    > IncrementalHydrostructure<'graph, 'walk, Graph, SafetyTracker>
{
    /// Compute the incremental hydrostructure of a walk.
    /// Sets the fingers to cover the whole walk.
    ///
    /// Panics if the given walk has less than two edges, since the hydrostructure is defined only for walks of at least two edges.
    pub fn compute(graph: &'graph Graph, walk: &'walk [Graph::EdgeIndex]) -> Self {
        debug_assert!(
            walk.len() >= 2,
            "The hydrostructure is defined only for walks of at least two edges."
        );
        let r_plus = compute_incremental_restricted_forward_edge_reachability(graph, walk);
        let r_minus = compute_incremental_restricted_backward_edge_reachability(graph, walk);
        let mut result = Self {
            left_finger: 0,
            right_finger: walk.len() - 1,
            rightmost_split: None,
            rightmost_join: None,
            safety_tracker: SafetyTracker::new_with_empty_subgraph(graph),
            r_plus,
            r_minus,
            walk,
        };
        result.reset_fingers();
        result
    }

    /// Compute the incremental hydrostructure of a walk.
    /// Sets the fingers to the given values.
    ///
    /// Panics if the given walk has less than two edges, since the hydrostructure is defined only for walks of at least two edges.
    /// Also panics if the fingers are set outside of the walk or the left finger is not left of the right finger.
    pub fn compute_and_set_fingers_to(
        graph: &'graph Graph,
        walk: &'walk [Graph::EdgeIndex],
        left_finger: usize,
        right_finger: usize,
    ) -> Self {
        debug_assert!(
            walk.len() >= 2,
            "The hydrostructure is defined only for walks of at least two edges."
        );
        debug_assert!(
            left_finger < right_finger,
            "Left finger must be smaller than the right finger."
        );
        debug_assert!(
            right_finger < walk.len(),
            "Thr right finger must be inside the walk."
        );
        let r_plus = compute_incremental_restricted_forward_edge_reachability(graph, walk);
        let r_minus = compute_incremental_restricted_backward_edge_reachability(graph, walk);
        let mut result = Self {
            left_finger,
            right_finger,
            rightmost_split: None,
            rightmost_join: None,
            safety_tracker: SafetyTracker::new_with_empty_subgraph(graph),
            r_plus,
            r_minus,
            walk,
        };
        result.reset_fingers();
        result
    }

    /// Compute the incremental hydrostructure of a walk.
    /// Sets the fingers to the first and second edge.
    ///
    /// Panics if the given walk has less than two edges, since the hydrostructure is defined only for walks of at least two edges.
    pub fn compute_and_set_fingers_left(
        graph: &'graph Graph,
        walk: &'walk [Graph::EdgeIndex],
    ) -> Self {
        debug_assert!(
            walk.len() >= 2,
            "The hydrostructure is defined only for walks of at least two edges."
        );
        let mut r_plus = compute_incremental_restricted_forward_edge_reachability(graph, walk);
        let mut r_minus = compute_incremental_restricted_backward_edge_reachability(graph, walk);
        r_plus.set_current_step(1);
        r_minus.set_current_step(walk.len() - 1);
        let mut result = Self {
            left_finger: 0,
            right_finger: 1,
            rightmost_split: if graph.is_split_edge(walk[1]) {
                Some(1)
            } else {
                None
            },
            rightmost_join: if graph.is_join_edge(walk[0]) {
                Some(0)
            } else {
                None
            },
            safety_tracker: SafetyTracker::new_with_empty_subgraph(graph),
            r_plus,
            r_minus,
            walk,
        };
        result.safety_tracker.reset(&result.r_plus, &result.r_minus);
        result
    }

    fn reset_rightmost_split_join(&mut self) {
        self.rightmost_split = self
            .walk
            .iter()
            .enumerate()
            .take(self.right_finger + 1)
            .skip(self.left_finger + 1)
            .rev()
            .filter(|(_, e)| self.r_plus.parent_graph().is_split_edge(**e))
            .map(|(n, _)| n)
            .next();
        self.rightmost_join = self
            .walk
            .iter()
            .enumerate()
            .take(self.right_finger)
            .skip(self.left_finger)
            .rev()
            .filter(|(_, e)| self.r_plus.parent_graph().is_join_edge(**e))
            .map(|(n, _)| n)
            .next();
    }

    /// When setting the fingers to arbitrary values, this computes the correct values for the rightmost split and join and the vapor_is_path_tracker.
    fn reset_fingers(&mut self) {
        self.r_minus
            .set_current_step(self.walk.len() - 1 - self.left_finger);
        self.r_plus.set_current_step(self.right_finger);

        self.reset_rightmost_split_join();
        self.safety_tracker.reset(&self.r_plus, &self.r_minus);
    }

    /// Set the left and right finger to the specified values.
    /// Panics if the given indices are not valid indices in the underlying walk of this incremental hydrostructure, or if the left finger is not left of the right finger.
    pub fn set_both_fingers(&mut self, left_finger: usize, right_finger: usize) {
        debug_assert!(
            left_finger < self.walk.len()
                && left_finger < right_finger
                && right_finger < self.walk.len()
        );
        self.left_finger = left_finger;
        self.right_finger = right_finger;
        self.reset_fingers();
    }

    /// Set the left finger to the specified value.
    /// Panics if the given index is not a valid index in the underlying walk of this incremental hydrostructure, or if it is not left of the right finger.
    pub fn set_left_finger(&mut self, left_finger: usize) {
        debug_assert!(left_finger < self.walk.len() && left_finger < self.right_finger);
        self.left_finger = left_finger;
        self.reset_fingers();
    }

    /// Returns the value of the left finger.
    pub fn left_finger(&self) -> usize {
        self.left_finger
    }

    /// Increments the left finger.
    /// Panics if the increment moves the left finger past the end of the walk.
    pub fn increment_left_finger(&mut self) {
        self.safety_tracker
            .remove_incremental_subgraph_step(&self.r_plus, &self.r_minus);

        self.left_finger += 1;
        debug_assert!(self.left_finger < self.walk.len() && self.left_finger < self.right_finger);
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

    /// Returns true if the left finger can be incremented without running into the right finger.
    pub fn can_increment_left_finger(&self) -> bool {
        self.left_finger + 1 < self.right_finger
    }

    /// Set the right finger to the specified value.
    /// Panics if the given index is not a valid index in the underlying walk of this incremental hydrostructure, or if it is not right of the left finger.
    pub fn set_right_finger(&mut self, right_finger: usize) {
        debug_assert!(right_finger < self.walk.len() && self.left_finger < right_finger);
        self.right_finger = right_finger;
        self.reset_fingers();
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
        debug_assert!(self.right_finger < self.walk.len());
        self.r_plus.set_current_step(self.right_finger);

        if self
            .r_plus
            .parent_graph()
            .is_split_edge(self.walk[self.right_finger])
        {
            self.rightmost_split = Some(self.right_finger);
        }

        self.safety_tracker
            .add_incremental_subgraph_step(&self.r_plus, &self.r_minus);
    }

    /// Returns true if the right finger can be incremented without running over the end of the walk.
    pub fn can_increment_right_finger(&self) -> bool {
        self.right_finger + 1 < self.walk.len()
    }

    /// Returns the subwalk from the left to the right finger (inclusive).
    pub fn current_walk(&self) -> &[Graph::EdgeIndex] {
        &self.walk[self.left_finger..self.right_finger + 1]
    }

    /// Returns true if the given node is in the path between left and right finger.
    fn is_node_in_trivial_r_plus_r_minus(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        self.walk
            .iter()
            .take(self.right_finger)
            .skip(self.left_finger)
            .map(|e| self.r_plus.parent_graph().edge_endpoints(*e).to_node)
            .any(|n| n == node)
    }

    /// Returns true if the given edge is any of the edges of the (inclusive) subwalk from left to right finger, except for the last one.
    fn is_edge_in_trivial_r_plus(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        self.walk
            .iter()
            .take(self.right_finger)
            .skip(self.left_finger)
            .any(|e| *e == edge)
    }

    /// Returns true if the given edge is any of the edges of the (inclusive) subwalk from left to right finger, except for the first one.
    fn is_edge_in_trivial_r_minus(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        self.walk
            .iter()
            .take(self.right_finger + 1)
            .skip(self.left_finger + 1)
            .any(|e| *e == edge)
    }

    /// Check if a node is in r_plus under the assumption that the current walk is bridge-like.
    fn is_node_r_plus_bridge_like(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        if self.rightmost_split.is_some() {
            self.r_plus.contains_node(node)
        } else {
            self.is_node_in_trivial_r_plus_r_minus(node)
        }
    }

    /// Check if a node is in r_minus under the assumption that the current walk is bridge-like.
    fn is_node_r_minus_bridge_like(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        if self.rightmost_join.is_some() {
            self.r_minus.contains_node(node)
        } else {
            self.is_node_in_trivial_r_plus_r_minus(node)
        }
    }

    /// Check if an edge is in r_plus under the assumption that the current walk is bridge-like.
    fn is_edge_r_plus_bridge_like(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        if self.rightmost_split.is_some() {
            DecoratingSubgraph::contains_edge(&self.r_plus, edge)
        } else {
            self.is_edge_in_trivial_r_plus(edge)
        }
    }

    /// Check if an edge is in r_minus under the assumption that the current walk is bridge-like.
    fn is_edge_r_minus_bridge_like(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        if self.rightmost_join.is_some() {
            DecoratingSubgraph::contains_edge(&self.r_minus, edge)
        } else {
            self.is_edge_in_trivial_r_minus(edge)
        }
    }

    /// Returns true if the current subwalk is safe according to the incremental safety tracker.
    pub fn is_safe(&self) -> bool {
        self.safety_tracker.is_safe(
            self.rightmost_split.is_none(),
            self.rightmost_join.is_none(),
        )
    }
}

impl<
        'graph,
        'walk,
        Graph: 'graph + StaticGraph,
        SafetyTracker: IncrementalSafetyTracker<'graph, Graph>,
    > Hydrostructure<Graph::NodeIndex, Graph::EdgeIndex>
    for IncrementalHydrostructure<'graph, 'walk, Graph, SafetyTracker>
{
    fn is_node_r_plus(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        if self.is_bridge_like() {
            self.is_node_r_plus_bridge_like(node)
        } else {
            true
        }
    }

    fn is_node_r_minus(&self, node: <Graph as GraphBase>::NodeIndex) -> bool {
        if self.is_bridge_like() {
            self.is_node_r_minus_bridge_like(node)
        } else {
            true
        }
    }

    fn is_edge_r_plus(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        if self.is_bridge_like() {
            self.is_edge_r_plus_bridge_like(edge)
        } else {
            true
        }
    }

    fn is_edge_r_minus(&self, edge: <Graph as GraphBase>::EdgeIndex) -> bool {
        if self.is_bridge_like() {
            self.is_edge_r_minus_bridge_like(edge)
        } else {
            true
        }
    }

    fn is_bridge_like(&self) -> bool {
        if SafetyTracker::does_safety_equal_bridge_like() {
            self.is_safe()
        } else {
            unimplemented!("The safety tracker used is not equivalent to a bridge-like check.")
        }
    }
}

/// An incremental safety tracker for the incremental hydrostructure.
pub trait IncrementalSafetyTracker<'a, Graph: GraphBase> {
    /// Creates a new instance keeping the hydrostructure components empty.
    fn new_with_empty_subgraph(graph: &'a Graph) -> Self;

    /// Remove all nodes and edges from the hydrostructure components.
    fn clear(&mut self);

    /// Reset the state of the safety tracker to the current state of the hydrostructure.
    fn reset(&mut self, r_plus: &IncrementalSubgraph<Graph>, r_minus: &IncrementalSubgraph<Graph>);

    /// Add the nodes and edges from the current step of r_plus.
    fn add_incremental_subgraph_step(
        &mut self,
        r_plus: &IncrementalSubgraph<Graph>,
        r_minus: &IncrementalSubgraph<Graph>,
    );

    /// Remove the nodes and edges from the current step of r_minus.
    fn remove_incremental_subgraph_step(
        &mut self,
        r_plus: &IncrementalSubgraph<Graph>,
        r_minus: &IncrementalSubgraph<Graph>,
    );

    /// Returns true if the safety tracker indicates that the current subwalk is safe.
    fn is_safe(&self, is_forward_univocal: bool, is_backward_univocal: bool) -> bool;

    /// May return true if the `is_safe` function returns true if and only if a subwalk is bridge-like.
    /// It is not required to return true for all types, but for types where it does return true
    /// the equivalence above must hold.
    fn does_safety_equal_bridge_like() -> bool;
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
    use crate::hydrostructure::incremental_hydrostructure::{BridgeLikeIncrementalHydrostructure};
    use crate::hydrostructure::static_hydrostructure::StaticHydrostructure;
    use crate::hydrostructure::Hydrostructure;
    use traitsequence::interface::Sequence;

    #[test]
    fn test_incremental_hydrostructure_path_2() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());

        let walk: Vec<_> = graph.create_edge_walk(&[e0, e1]);
        let incremental_hydrostructure =
            BridgeLikeIncrementalHydrostructure::compute(&graph, &walk);
        debug_assert!(incremental_hydrostructure.is_node_river(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e1));
        debug_assert!(incremental_hydrostructure.is_node_river(n2));
    }

    #[test]
    fn test_incremental_hydrostructure_path_3() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n2, n3, ());

        let walk: Vec<_> = graph.create_edge_walk(&[e0, e1, e2]);
        let mut incremental_hydrostructure =
            BridgeLikeIncrementalHydrostructure::compute(&graph, &walk);
        debug_assert!(incremental_hydrostructure.is_node_river(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_vapor(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));
        debug_assert!(incremental_hydrostructure.is_node_river(n3));

        incremental_hydrostructure.set_both_fingers(0, 1);
        debug_assert!(incremental_hydrostructure.is_node_river(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e1));
        debug_assert!(incremental_hydrostructure.is_node_river(n2));
        debug_assert!(incremental_hydrostructure.is_edge_river(e2));
        debug_assert!(incremental_hydrostructure.is_node_river(n3));

        incremental_hydrostructure.increment_right_finger();
        debug_assert!(incremental_hydrostructure.is_node_river(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_vapor(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));
        debug_assert!(incremental_hydrostructure.is_node_river(n3));

        incremental_hydrostructure.increment_left_finger();
        debug_assert!(incremental_hydrostructure.is_node_river(n0));
        debug_assert!(incremental_hydrostructure.is_edge_river(e0));
        debug_assert!(incremental_hydrostructure.is_node_river(n1));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));
        debug_assert!(incremental_hydrostructure.is_node_river(n3));
    }

    #[test]
    fn test_incremental_hydrostructure_cycle_start_with_loop() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let e0 = graph.add_edge(n0, n0, ());
        let e1 = graph.add_edge(n0, n1, ());
        let e2 = graph.add_edge(n1, n2, ());
        let e3 = graph.add_edge(n2, n0, ());

        let walk: Vec<_> = graph.create_edge_walk(&[e0, e1, e2]);
        let mut incremental_hydrostructure =
            BridgeLikeIncrementalHydrostructure::compute(&graph, &walk);
        debug_assert!(incremental_hydrostructure.is_node_vapor(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_edge_vapor(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));
        debug_assert!(incremental_hydrostructure.is_node_cloud(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e3));

        incremental_hydrostructure.set_both_fingers(0, 1);
        debug_assert!(incremental_hydrostructure.is_node_vapor(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e1));
        debug_assert!(incremental_hydrostructure.is_node_cloud(n1));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));
        debug_assert!(incremental_hydrostructure.is_node_cloud(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e3));

        incremental_hydrostructure.increment_right_finger();
        debug_assert!(incremental_hydrostructure.is_node_vapor(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_edge_vapor(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));
        debug_assert!(incremental_hydrostructure.is_node_cloud(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e3));

        incremental_hydrostructure.increment_left_finger();
        debug_assert!(incremental_hydrostructure.is_node_river(n0));
        debug_assert!(incremental_hydrostructure.is_edge_river(e0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));
        debug_assert!(incremental_hydrostructure.is_node_river(n2));
        debug_assert!(incremental_hydrostructure.is_edge_river(e3));
    }

    #[test]
    fn test_incremental_hydrostructure_cycle_end_with_loop() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let e0 = graph.add_edge(n0, n1, ());
        let e1 = graph.add_edge(n1, n2, ());
        let e2 = graph.add_edge(n2, n2, ());
        let e3 = graph.add_edge(n2, n0, ());

        let walk: Vec<_> = graph.create_edge_walk(&[e0, e1, e2]);
        let mut incremental_hydrostructure =
            BridgeLikeIncrementalHydrostructure::compute(&graph, &walk);
        debug_assert!(incremental_hydrostructure.is_edge_sea(e3));
        debug_assert!(incremental_hydrostructure.is_node_sea(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_vapor(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));

        incremental_hydrostructure.set_both_fingers(0, 1);
        debug_assert!(incremental_hydrostructure.is_edge_river(e3));
        debug_assert!(incremental_hydrostructure.is_node_river(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e1));
        debug_assert!(incremental_hydrostructure.is_node_river(n2));
        debug_assert!(incremental_hydrostructure.is_edge_river(e2));

        incremental_hydrostructure.increment_right_finger();
        debug_assert!(incremental_hydrostructure.is_edge_sea(e3));
        debug_assert!(incremental_hydrostructure.is_node_sea(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n1));
        debug_assert!(incremental_hydrostructure.is_edge_vapor(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));

        incremental_hydrostructure.increment_left_finger();
        debug_assert!(incremental_hydrostructure.is_edge_sea(e3));
        debug_assert!(incremental_hydrostructure.is_node_sea(n0));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e0));
        debug_assert!(incremental_hydrostructure.is_node_sea(n1));
        debug_assert!(incremental_hydrostructure.is_edge_sea(e1));
        debug_assert!(incremental_hydrostructure.is_node_vapor(n2));
        debug_assert!(incremental_hydrostructure.is_edge_cloud(e2));
    }

    #[test]
    fn test_incremental_hydrostructure_macrotig_bridge_like_only() {
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

        for macrotig in maximal_macrotigs.iter() {
            let mut incremental_hydrostructure =
                BridgeLikeIncrementalHydrostructure::compute(&graph, macrotig);
            for subwalk_len in 2..macrotig.len() {
                //println!("Setting initial fingers for subwalk len {}", subwalk_len);
                incremental_hydrostructure.set_both_fingers(0, subwalk_len - 1);
                for offset in 0..(macrotig.len() - subwalk_len) {
                    //let subwalk: VecEdgeWalk<_> =
                    //    graph.create_edge_walk(&macrotig[offset..offset + subwalk_len]);
                    let subwalk = Vec::from(&macrotig[offset..offset + subwalk_len]);
                    //println!("{:?}", subwalk);
                    /*println!(
                        "has rightmost split/join: {:?}/{:?}; left/right finger: {}/{}",
                        incremental_hydrostructure.rightmost_split,
                        incremental_hydrostructure.rightmost_join,
                        incremental_hydrostructure.left_finger,
                        incremental_hydrostructure.right_finger
                    );*/
                    let static_hydrostructure =
                        StaticHydrostructure::compute_with_bitvector_subgraph(
                            &graph,
                            subwalk.clone(),
                        );

                    if static_hydrostructure.is_bridge_like() {
                        for node in graph.node_indices() {
                            debug_assert_eq!(
                                incremental_hydrostructure.is_node_r_plus(node),
                                static_hydrostructure.is_node_r_plus(node),
                                "Node {:?} r_plus incremental/static: {}/{}",
                                node,
                                incremental_hydrostructure.is_node_r_plus(node),
                                static_hydrostructure.is_node_r_plus(node)
                            );
                            debug_assert_eq!(
                                incremental_hydrostructure.is_node_r_minus(node),
                                static_hydrostructure.is_node_r_minus(node),
                                "Node {:?} r_minus incremental/static: {}/{}",
                                node,
                                incremental_hydrostructure.is_node_r_minus(node),
                                static_hydrostructure.is_node_r_minus(node)
                            );
                        }
                        for edge in graph.edge_indices() {
                            debug_assert_eq!(
                                incremental_hydrostructure.is_edge_r_plus(edge),
                                static_hydrostructure.is_edge_r_plus(edge),
                                "Edge {:?} r_plus incremental/static: {}/{}",
                                edge,
                                incremental_hydrostructure.is_edge_r_plus(edge),
                                static_hydrostructure.is_edge_r_plus(edge)
                            );
                            debug_assert_eq!(
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
                        /*println!(
                            "Before incrementing fingers: {:?}",
                            incremental_hydrostructure.vapor_is_path_tracker
                        );*/
                        incremental_hydrostructure.increment_right_finger();
                        /*println!(
                            "Before incrementing left finger: {:?}",
                            incremental_hydrostructure.vapor_is_path_tracker
                        );*/
                        incremental_hydrostructure.increment_left_finger();
                    }
                }
            }
        }
    }

    #[test]
    fn test_incremental_hydrostructure_macrotig_bridge_like_and_avertible() {
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

        for macrotig in maximal_macrotigs.iter() {
            let mut incremental_hydrostructure =
                BridgeLikeIncrementalHydrostructure::compute(&graph, macrotig);
            for subwalk_len in 2..macrotig.len() {
                //println!("Setting initial left finger for subwalk len {}", subwalk_len);
                incremental_hydrostructure.set_left_finger(0);
                //println!("Setting initial right finger for subwalk len {}", subwalk_len);
                incremental_hydrostructure.set_right_finger(subwalk_len - 1);
                for offset in 0..(macrotig.len() - subwalk_len) {
                    //let subwalk: VecEdgeWalk<_> =
                    //    graph.create_edge_walk(&macrotig[offset..offset + subwalk_len]);
                    let subwalk = Vec::from(&macrotig[offset..offset + subwalk_len]);
                    //println!("{:?}", subwalk);
                    //println!("has rightmost split/join: {:?}/{:?}; left/right finger: {}/{}", incremental_hydrostructure.rightmost_split, incremental_hydrostructure.rightmost_join, incremental_hydrostructure.left_finger, incremental_hydrostructure.right_finger);
                    let static_hydrostructure =
                        StaticHydrostructure::compute_with_bitvector_subgraph(
                            &graph,
                            subwalk.clone(),
                        );

                    for node in graph.node_indices() {
                        debug_assert_eq!(
                            incremental_hydrostructure.is_node_r_plus(node),
                            static_hydrostructure.is_node_r_plus(node),
                            "Node {:?} r_plus incremental/static: {}/{}",
                            node,
                            incremental_hydrostructure.is_node_r_plus(node),
                            static_hydrostructure.is_node_r_plus(node)
                        );
                        debug_assert_eq!(
                            incremental_hydrostructure.is_node_r_minus(node),
                            static_hydrostructure.is_node_r_minus(node),
                            "Node {:?} r_minus incremental/static: {}/{}",
                            node,
                            incremental_hydrostructure.is_node_r_minus(node),
                            static_hydrostructure.is_node_r_minus(node)
                        );
                    }
                    for edge in graph.edge_indices() {
                        debug_assert_eq!(
                            incremental_hydrostructure.is_edge_r_plus(edge),
                            static_hydrostructure.is_edge_r_plus(edge),
                            "Edge {:?} r_plus incremental/static: {}/{}",
                            edge,
                            incremental_hydrostructure.is_edge_r_plus(edge),
                            static_hydrostructure.is_edge_r_plus(edge)
                        );
                        debug_assert_eq!(
                            incremental_hydrostructure.is_edge_r_minus(edge),
                            static_hydrostructure.is_edge_r_minus(edge),
                            "Edge {:?} r_minus incremental/static: {}/{}",
                            edge,
                            incremental_hydrostructure.is_edge_r_minus(edge),
                            static_hydrostructure.is_edge_r_minus(edge)
                        );
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
