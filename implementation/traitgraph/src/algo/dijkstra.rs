use crate::index::GraphIndex;
use crate::interface::{GraphBase, StaticGraph};
use std::collections::{BinaryHeap, HashMap, HashSet};

/// A Dijkstra implementation with a set of common optimisations.
pub type DefaultDijkstra<'a, Graph> = Dijkstra<'a, Graph, EpochNodeWeightArray<usize>>;

/// A weight-type usable in Dijkstra's algorithm.
pub trait Weight {
    /// The infinity value of this type.
    fn infinity() -> Self;
}

impl Weight for usize {
    #[inline]
    fn infinity() -> Self {
        Self::max_value()
    }
}

/// Edge data that has a weight usable for shortest path computation.
pub trait WeightedEdgeData {
    /// The weight of the edge.
    fn weight(&self) -> usize;
}

impl WeightedEdgeData for usize {
    #[inline]
    fn weight(&self) -> usize {
        *self
    }
}

/// An array to store minimal node weights for Dijkstra's algorithm.
pub trait NodeWeightArray<WeightType> {
    /// Create a new NodeWeightArray of given size.
    fn new(size: usize) -> Self;

    /// Returns the current weight of the given node index.
    fn get(&self, node_index: usize) -> WeightType;

    /// Returns the current weight of the given node index as mutable reference.
    fn get_mut(&mut self, node_index: usize) -> &mut WeightType;

    /// Sets the current weight of the given node index.
    fn set(&mut self, node_index: usize, weight: WeightType);

    /// Resets the weights of all node indices to infinity
    fn clear(&mut self);
}

impl<WeightType: Weight + Copy> NodeWeightArray<WeightType> for Vec<WeightType> {
    fn new(size: usize) -> Self {
        vec![WeightType::infinity(); size]
    }

    #[inline]
    fn get(&self, node_index: usize) -> WeightType {
        self[node_index]
    }

    #[inline]
    fn get_mut(&mut self, node_index: usize) -> &mut WeightType {
        &mut self[node_index]
    }

    #[inline]
    fn set(&mut self, node_index: usize, weight: WeightType) {
        self[node_index] = weight;
    }

    fn clear(&mut self) {
        for entry in self.iter_mut() {
            *entry = WeightType::infinity();
        }
    }
}

/// An epoch counter array.
/// This can be used to check if an index is current by comparing its entry in the epoch array to the current epoch.
/// To unmark all values, the current epoch can be increased in O(1). Only overflows have to be handled by resetting all epoch counters.
pub struct EpochArray {
    epochs: Vec<u32>,
    current_epoch: u32,
}

impl EpochArray {
    /// Create a new epoch array of given length where all values are outdated.
    pub fn new(len: usize) -> Self {
        Self {
            epochs: vec![0; len],
            current_epoch: 1,
        }
    }

    /// Outdate all indices.
    pub fn clear(&mut self) {
        if self.current_epoch == u32::max_value() {
            for epoch in self.epochs.iter_mut() {
                *epoch = 0;
            }
            self.current_epoch = 1;
        } else {
            self.current_epoch += 1;
        }
    }

    /// Set the given index as current.
    #[inline]
    pub fn update(&mut self, index: usize) {
        self.epochs[index] = self.current_epoch;
    }

    /// Returns true if the given index is current, and false otherwise.
    #[inline]
    pub fn get(&self, index: usize) -> bool {
        self.epochs[index] == self.current_epoch
    }

    /// Updates the given index and returns true if the given index was current before, and false otherwise.
    #[inline]
    pub fn get_and_update(&mut self, index: usize) -> bool {
        if self.epochs[index] == self.current_epoch {
            true
        } else {
            self.epochs[index] = self.current_epoch;
            false
        }
    }
}

/// An epoched node weight array that can be cleared in O(1) most of the times.
/// Only if the epoch in the epoch array overflows, clearing takes linear time.
pub struct EpochNodeWeightArray<WeightType> {
    weights: Vec<WeightType>,
    epochs: EpochArray,
}

impl<WeightType: Weight> EpochNodeWeightArray<WeightType> {
    #[inline]
    fn make_current(&mut self, node_index: usize) {
        if !self.epochs.get_and_update(node_index) {
            self.weights[node_index] = WeightType::infinity();
        }
    }
}

impl<WeightType: Weight + Copy> NodeWeightArray<WeightType> for EpochNodeWeightArray<WeightType> {
    fn new(len: usize) -> Self {
        Self {
            weights: vec![WeightType::infinity(); len],
            epochs: EpochArray::new(len),
        }
    }

    #[inline]
    fn get(&self, node_index: usize) -> WeightType {
        if self.epochs.get(node_index) {
            self.weights[node_index]
        } else {
            WeightType::infinity()
        }
    }

    #[inline]
    fn get_mut(&mut self, node_index: usize) -> &mut WeightType {
        self.make_current(node_index);
        &mut self.weights[node_index]
    }

    #[inline]
    fn set(&mut self, node_index: usize, weight: WeightType) {
        self.weights[node_index] = weight;
        self.epochs.update(node_index);
    }

    fn clear(&mut self) {
        self.epochs.clear();
    }
}

/// Datastructure for Dijkstra's shortest path algorithm.
pub struct Dijkstra<'a, Graph: GraphBase, NodeWeights> {
    queue: BinaryHeap<std::cmp::Reverse<(usize, Graph::NodeIndex)>>,
    // back_pointers: Vec<Graph::OptionalNodeIndex>,
    node_weights: NodeWeights,
    target_nodes: EpochArray,
    graph: &'a Graph,
}

impl<
        'a,
        EdgeData: WeightedEdgeData,
        Graph: StaticGraph<EdgeData = EdgeData>,
        NodeWeights: NodeWeightArray<usize>,
    > Dijkstra<'a, Graph, NodeWeights>
{
    /// Create the datastructures for the given graph.
    pub fn new(graph: &'a Graph) -> Self {
        Self {
            queue: BinaryHeap::new(),
            // back_pointers: vec![Default::default(); graph.node_count()],
            node_weights: NodeWeights::new(graph.node_count()),
            target_nodes: EpochArray::new(graph.node_count()),
            graph,
        }
    }

    /// Compute the shortest paths from source to all targets, with given maximum weight.
    #[inline(never)]
    pub fn shortest_path_lens(
        &mut self,
        source: Graph::NodeIndex,
        targets: &HashSet<Graph::NodeIndex>,
        max_weight: usize,
    ) -> HashMap<Graph::NodeIndex, usize> {
        //println!("Shortest path lens of {}", source.as_usize());
        let mut shortest_path_lens = HashMap::new();

        self.queue.push(std::cmp::Reverse((0, source)));
        //self.back_pointers[source.as_usize()] = source.into();
        self.node_weights.set(source.as_usize(), 0);

        for target in targets {
            self.target_nodes.update(target.as_usize());
        }

        while let Some(std::cmp::Reverse((weight, node_index))) = self.queue.pop() {
            //println!("Finalising node {}", node_index.as_usize());
            // Check if the node was already processed
            let actual_weight = self.node_weights.get(node_index.as_usize());
            if actual_weight < weight {
                continue;
            }
            assert_eq!(actual_weight, weight);

            // Check if we already found all paths
            if shortest_path_lens.len() == targets.len() {
                break;
            }

            // Check if we are still below max_weight
            if weight > max_weight {
                break;
            }

            // Check if we found a target
            if self.target_nodes.get(node_index.as_usize()) {
                shortest_path_lens.insert(node_index, weight);
            }
            /*if targets.contains(&node_index) {
                shortest_path_lens.insert(node_index, weight);
            }*/

            // Relax neighbors
            for out_neighbor in self.graph.out_neighbors(node_index) {
                let new_neighbor_weight =
                    weight + self.graph.edge_data(out_neighbor.edge_id).weight();
                let neighbor_weight = self.node_weights.get_mut(out_neighbor.node_id.as_usize());
                if new_neighbor_weight < *neighbor_weight {
                    *neighbor_weight = new_neighbor_weight;
                    self.queue.push(std::cmp::Reverse((
                        new_neighbor_weight,
                        out_neighbor.node_id,
                    )));
                    //self.back_pointers[out_neighbor.node_id.as_usize()] = node_index.into();
                }
            }
        }

        self.queue.clear();
        /*for back_pointer in &mut self.back_pointers {
            *back_pointer = Default::default();
        }*/
        self.node_weights.clear();
        self.target_nodes.clear();
        shortest_path_lens
    }
}

#[cfg(test)]
mod tests {
    use crate::algo::dijkstra::DefaultDijkstra;
    use crate::implementation::petgraph_impl;
    use crate::interface::MutableGraphContainer;

    #[test]
    fn test_dijkstra_simple() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        graph.add_edge(n1, n2, 2);
        graph.add_edge(n2, n3, 2);
        graph.add_edge(n1, n3, 5);

        let mut dijkstra = DefaultDijkstra::new(&graph);
        let shortest_path_lens =
            dijkstra.shortest_path_lens(n1, &[n3].iter().copied().collect(), 6);
        assert_eq!(shortest_path_lens, [(n3, 4)].iter().copied().collect());

        let shortest_path_lens =
            dijkstra.shortest_path_lens(n1, &[n3].iter().copied().collect(), 6);
        assert_eq!(shortest_path_lens, [(n3, 4)].iter().copied().collect());

        let shortest_path_lens =
            dijkstra.shortest_path_lens(n2, &[n3].iter().copied().collect(), 6);
        assert_eq!(shortest_path_lens, [(n3, 2)].iter().copied().collect());

        let shortest_path_lens =
            dijkstra.shortest_path_lens(n3, &[n3].iter().copied().collect(), 6);
        assert_eq!(shortest_path_lens, [(n3, 0)].iter().copied().collect());

        let shortest_path_lens =
            dijkstra.shortest_path_lens(n3, &[n2].iter().copied().collect(), 6);
        assert_eq!(shortest_path_lens, [].iter().copied().collect());
    }

    #[test]
    fn test_dijkstra_cycle() {
        let mut graph = petgraph_impl::new();
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        graph.add_edge(n1, n2, 2);
        graph.add_edge(n2, n3, 2);
        graph.add_edge(n3, n1, 5);

        let mut dijkstra = DefaultDijkstra::new(&graph);
        let shortest_path_lens =
            dijkstra.shortest_path_lens(n1, &[n3].iter().copied().collect(), 6);
        assert_eq!(shortest_path_lens, [(n3, 4)].iter().copied().collect());
    }
}
