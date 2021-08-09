use crate::index::GraphIndex;
use crate::interface::{GraphBase, StaticGraph};
use std::collections::BinaryHeap;
use std::marker::PhantomData;

/// A Dijkstra implementation with a set of common optimisations.
pub type DefaultDijkstra<Graph> = Dijkstra<Graph, EpochNodeWeightArray<usize>>;
//pub type DefaultDijkstra<'a, Graph> = Dijkstra<'a, Graph, Vec<usize>>;

/// A weight-type usable in Dijkstra's algorithm.
pub trait Weight {
    /// The infinity value of this type.
    fn infinity() -> Self;
}

impl Weight for usize {
    #[inline]
    fn infinity() -> Self {
        Self::MAX
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
        unsafe {
            *self.epochs.get_unchecked_mut(index) = self.current_epoch;
        }
        //self.epochs[index] = self.current_epoch;
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
pub struct Dijkstra<Graph: GraphBase, NodeWeights> {
    queue: BinaryHeap<std::cmp::Reverse<(usize, Graph::NodeIndex)>>,
    // back_pointers: Vec<Graph::OptionalNodeIndex>,
    node_weights: NodeWeights,
    graph: PhantomData<Graph>,
}

impl<
        EdgeData: WeightedEdgeData,
        Graph: StaticGraph<EdgeData = EdgeData>,
        NodeWeights: NodeWeightArray<usize>,
    > Dijkstra<Graph, NodeWeights>
{
    /// Create the datastructures for the given graph.
    pub fn new(graph: &Graph) -> Self {
        Self {
            queue: BinaryHeap::new(),
            // back_pointers: vec![Default::default(); graph.node_count()],
            node_weights: NodeWeights::new(graph.node_count()),
            graph: Default::default(),
        }
    }

    /// Compute the shortest paths from source to all targets, with given maximum weight.
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    pub fn shortest_path_lens<TargetMap: std::ops::Index<usize, Output = bool>>(
        &mut self,
        graph: &Graph,
        source: Graph::NodeIndex,
        targets: &TargetMap,
        target_amount: usize,
        max_weight: usize,
        forbid_source_target: bool,
        distances: &mut Vec<(Graph::NodeIndex, usize)>,
    ) {
        //println!("Shortest path lens of {}", source.as_usize());
        self.queue.push(std::cmp::Reverse((0, source)));
        //self.back_pointers[source.as_usize()] = source.into();
        self.node_weights.set(source.as_usize(), 0);
        distances.clear();

        //let mut iterations = 0;
        //let mut unnecessary_iterations = 0;
        //let max_iterations = self.graph.node_count();
        while let Some(std::cmp::Reverse((weight, node_index))) = self.queue.pop() {
            //iterations += 1;
            //println!("Finalising node {}", node_index.as_usize());
            // Check if the node was already processed
            let actual_weight = self.node_weights.get(node_index.as_usize());
            if actual_weight < weight {
                //unnecessary_iterations += 1;
                continue;
            }
            assert_eq!(actual_weight, weight);

            // Check if we are still lower than or equal to max_weight
            if weight > max_weight {
                //println!("Aborting early by max_weight after {}/{} iterations of which {} are unnecessary", iterations, max_iterations, unnecessary_iterations);
                break;
            }

            // Check if we found a target
            if targets[node_index.as_usize()] && (!forbid_source_target || node_index != source) {
                distances.push((node_index, weight));

                // Check if we already found all paths
                if distances.len() == target_amount {
                    //println!("Aborting early after finding all targets");
                    break;
                }
            }

            // Relax neighbors
            for out_neighbor in graph.out_neighbors(node_index) {
                let new_neighbor_weight = weight + graph.edge_data(out_neighbor.edge_id).weight();
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
        let mut distances = Vec::new();
        let mut targets = vec![false, false, true];
        dijkstra.shortest_path_lens(&graph, n1, &targets, 1, 6, false, &mut distances);
        assert_eq!(distances, vec![(n3, 4)]);

        dijkstra.shortest_path_lens(&graph, n1, &targets, 1, 6, false, &mut distances);
        assert_eq!(distances, vec![(n3, 4)]);

        dijkstra.shortest_path_lens(&graph, n2, &targets, 1, 6, false, &mut distances);
        assert_eq!(distances, vec![(n3, 2)]);

        dijkstra.shortest_path_lens(&graph, n3, &targets, 1, 6, false, &mut distances);
        assert_eq!(distances, vec![(n3, 0)]);

        targets = vec![false, true, false];
        dijkstra.shortest_path_lens(&graph, n3, &targets, 1, 6, false, &mut distances);
        assert_eq!(distances, vec![]);
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
        let mut distances = Vec::new();
        let targets = vec![false, false, true];
        dijkstra.shortest_path_lens(&graph, n1, &targets, 1, 6, false, &mut distances);
        assert_eq!(distances, vec![(n3, 4)]);
    }
}
