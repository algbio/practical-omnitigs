use crate::hydrostructure::incremental_hydrostructure::IncrementalSafetyTracker;
use traitgraph::implementation::incremental_subgraph::IncrementalSubgraph;
use traitgraph::interface::GraphBase;

/// An incremental safety tracker that tracks safety using the conjunction of two other safety trackers.
pub struct ConjunctiveSafetyTracker<SafetyTracker1, SafetyTracker2> {
    safety_tracker_1: SafetyTracker1,
    safety_tracker_2: SafetyTracker2,
}

impl<
        'a,
        Graph: GraphBase,
        SafetyTracker1: IncrementalSafetyTracker<'a, Graph>,
        SafetyTracker2: IncrementalSafetyTracker<'a, Graph>,
    > IncrementalSafetyTracker<'a, Graph>
    for ConjunctiveSafetyTracker<SafetyTracker1, SafetyTracker2>
{
    fn new_with_empty_subgraph(graph: &'a Graph) -> Self {
        Self {
            safety_tracker_1: SafetyTracker1::new_with_empty_subgraph(graph),
            safety_tracker_2: SafetyTracker2::new_with_empty_subgraph(graph),
        }
    }

    fn clear(&mut self) {
        self.safety_tracker_1.clear();
        self.safety_tracker_2.clear();
    }

    fn reset(&mut self, r_plus: &IncrementalSubgraph<Graph>, r_minus: &IncrementalSubgraph<Graph>) {
        self.safety_tracker_1.reset(r_plus, r_minus);
        self.safety_tracker_2.reset(r_plus, r_minus);
    }

    fn add_incremental_subgraph_step(
        &mut self,
        r_plus: &IncrementalSubgraph<Graph>,
        r_minus: &IncrementalSubgraph<Graph>,
    ) {
        self.safety_tracker_1
            .add_incremental_subgraph_step(r_plus, r_minus);
        self.safety_tracker_2
            .add_incremental_subgraph_step(r_plus, r_minus);
    }

    fn remove_incremental_subgraph_step(
        &mut self,
        r_plus: &IncrementalSubgraph<Graph>,
        r_minus: &IncrementalSubgraph<Graph>,
    ) {
        self.safety_tracker_1
            .remove_incremental_subgraph_step(r_plus, r_minus);
        self.safety_tracker_2
            .remove_incremental_subgraph_step(r_plus, r_minus);
    }

    fn is_safe(&self, is_forward_univocal: bool, is_backward_univocal: bool) -> bool {
        self.safety_tracker_1
            .is_safe(is_forward_univocal, is_backward_univocal)
            && self
                .safety_tracker_2
                .is_safe(is_forward_univocal, is_backward_univocal)
    }

    fn does_safety_equal_bridge_like() -> bool {
        SafetyTracker1::does_safety_equal_bridge_like()
            && SafetyTracker2::does_safety_equal_bridge_like()
    }
}
