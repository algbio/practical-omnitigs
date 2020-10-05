use crate::hydrostructure::incremental_hydrostructure::IncrementalSafetyTracker;
use traitgraph::implementation::incremental_subgraph::IncrementalSubgraph;
use traitgraph::interface::subgraph::DecoratingSubgraph;
use traitgraph::interface::ImmutableGraphContainer;

/// A type that keeps counts of the nodes in the different hydrostructure components to dynamically determine if they contain nodes.
/// It reports if a bridge-like walk is safe in the 1-circular node-centric model.
pub struct NodeCentricComponentTracker {
    sea_node_count: usize,
    vapor_node_count: usize,
    cloud_node_count: usize,
    river_node_count: usize,
}

impl<'a, Graph: ImmutableGraphContainer> IncrementalSafetyTracker<'a, Graph>
    for NodeCentricComponentTracker
{
    fn new_with_empty_subgraph(_graph: &'a Graph) -> Self {
        Self {
            sea_node_count: 0,
            vapor_node_count: 0,
            cloud_node_count: 0,
            river_node_count: 0,
        }
    }

    fn clear(&mut self) {
        self.sea_node_count = 0;
        self.vapor_node_count = 0;
        self.cloud_node_count = 0;
        self.river_node_count = 0;
    }

    fn reset(&mut self, r_plus: &IncrementalSubgraph<Graph>, r_minus: &IncrementalSubgraph<Graph>) {
        IncrementalSafetyTracker::<'a, Graph>::clear(self);
        for node in r_plus.parent_graph().node_indices() {
            match (r_plus.contains_node(node), r_minus.contains_node(node)) {
                (true, true) => {
                    self.vapor_node_count = self
                        .vapor_node_count
                        .checked_add(1)
                        .expect("Overflow in node count")
                }
                (true, false) => {
                    self.sea_node_count = self
                        .sea_node_count
                        .checked_add(1)
                        .expect("Overflow in node count")
                }
                (false, true) => {
                    self.cloud_node_count = self
                        .cloud_node_count
                        .checked_add(1)
                        .expect("Overflow in node count")
                }
                (false, false) => {
                    self.river_node_count = self
                        .river_node_count
                        .checked_add(1)
                        .expect("Overflow in node count")
                }
            }
        }
    }

    fn add_incremental_subgraph_step(
        &mut self,
        r_plus: &IncrementalSubgraph<Graph>,
        r_minus: &IncrementalSubgraph<Graph>,
    ) {
        for node in r_plus.new_nodes() {
            if r_minus.contains_node(*node) {
                self.cloud_node_count = self
                    .cloud_node_count
                    .checked_sub(1)
                    .expect("Overflow in node count");
                self.vapor_node_count = self
                    .vapor_node_count
                    .checked_add(1)
                    .expect("Overflow in node count");
            } else {
                self.river_node_count = self
                    .river_node_count
                    .checked_sub(1)
                    .expect("Overflow in node count");
                self.sea_node_count = self
                    .sea_node_count
                    .checked_add(1)
                    .expect("Overflow in node count");
            }
        }
    }

    fn remove_incremental_subgraph_step(
        &mut self,
        r_plus: &IncrementalSubgraph<Graph>,
        r_minus: &IncrementalSubgraph<Graph>,
    ) {
        for node in r_minus.new_nodes() {
            if r_plus.contains_node(*node) {
                self.vapor_node_count = self
                    .vapor_node_count
                    .checked_sub(1)
                    .expect("Overflow in node count");
                self.sea_node_count = self
                    .sea_node_count
                    .checked_add(1)
                    .expect("Overflow in node count");
            } else {
                self.cloud_node_count = self
                    .cloud_node_count
                    .checked_sub(1)
                    .expect("Overflow in node count");
                self.river_node_count = self
                    .river_node_count
                    .checked_add(1)
                    .expect("Overflow in node count");
            }
        }
    }

    fn is_safe(&self, is_forward_univocal: bool, is_backward_univocal: bool) -> bool {
        // We assume that the walk is bridge-like.
        if is_forward_univocal && is_backward_univocal {
            // If the walk is biunivocal then the head of its end is a node in the river.
            true
        } else if is_forward_univocal {
            // If the walk is forward-univocal, then all nodes in the sea fall into the river temporarily.
            // Also, the sea is then empty, so the sea-cloud condition does not hold.
            self.river_node_count + self.sea_node_count > 0
        } else if is_backward_univocal {
            // If the walk is backward-univocal, then all nodes in the cloud fall into the river temporarily.
            // Also, the cloud is then empty, so the sea-cloud condition does not hold.
            self.river_node_count + self.cloud_node_count > 0
        } else {
            self.river_node_count > 0 || (self.cloud_node_count > 0 && self.sea_node_count > 0)
        }
    }

    fn does_safety_equal_bridge_like() -> bool {
        false
    }
}
