use crate::index::GraphIndex;
use crate::interface::{GraphBase, StaticGraph};
use std::collections::{BinaryHeap, HashMap, HashSet};

/// Edge data that has a weight usable for shortest path computation.
pub trait WeightedEdgeData {
    /// The weight of the edge.
    fn weight(&self) -> usize;
}

impl WeightedEdgeData for usize {
    fn weight(&self) -> usize {
        *self
    }
}

/// Datastructure for Dijkstra's shortest path algorithm.
pub struct Dijkstra<'a, Graph: GraphBase> {
    queue: BinaryHeap<std::cmp::Reverse<(usize, Graph::NodeIndex)>>,
    // back_pointers: Vec<Graph::OptionalNodeIndex>,
    node_weights: Vec<usize>,
    graph: &'a Graph,
}

impl<'a, EdgeData: WeightedEdgeData, Graph: StaticGraph<EdgeData = EdgeData>> Dijkstra<'a, Graph> {
    /// Create the datastructures for the given graph.
    pub fn new(graph: &'a Graph) -> Self {
        Self {
            queue: BinaryHeap::new(),
            // back_pointers: vec![Default::default(); graph.node_count()],
            node_weights: vec![usize::max_value(); graph.node_count()],
            graph,
        }
    }

    /// Compute the shortest paths from source to all targets, with given maximum weight.
    pub fn shortest_path_lens(
        &mut self,
        source: Graph::NodeIndex,
        targets: &HashSet<Graph::NodeIndex>,
        max_weight: usize,
    ) -> HashMap<Graph::NodeIndex, usize> {
        let mut shortest_path_lens = HashMap::new();

        self.queue.push(std::cmp::Reverse((0, source)));
        //self.back_pointers[source.as_usize()] = source.into();
        self.node_weights[source.as_usize()] = 0;

        while let Some(std::cmp::Reverse((weight, node_index))) = self.queue.pop() {
            // Check if the node was already processed
            if self.node_weights[node_index.as_usize()] < weight {
                continue;
            }
            assert_eq!(self.node_weights[node_index.as_usize()], weight);

            // Check if we already found all paths
            if shortest_path_lens.len() == targets.len() {
                break;
            }

            // Check if we are still below max_weight
            if weight > max_weight {
                break;
            }

            // Check if we found a target
            if targets.contains(&node_index) {
                shortest_path_lens.insert(node_index, weight);
            }

            // Relax neighbors
            for out_neighbor in self.graph.out_neighbors(node_index) {
                let new_neighbor_weight =
                    weight + self.graph.edge_data(out_neighbor.edge_id).weight();
                if new_neighbor_weight < self.node_weights[out_neighbor.node_id.as_usize()] {
                    self.node_weights[out_neighbor.node_id.as_usize()] = new_neighbor_weight;
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
        for node_weight in &mut self.node_weights {
            *node_weight = usize::max_value();
        }
        shortest_path_lens
    }
}

#[cfg(test)]
mod tests {
    use crate::algo::dijkstra::Dijkstra;
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

        let mut dijkstra = Dijkstra::new(&graph);
        let shortest_path_lens =
            dijkstra.shortest_path_lens(n1, &[n3].iter().copied().collect(), 6);
        assert_eq!(shortest_path_lens, [(n3, 4)].iter().copied().collect());
    }
}
