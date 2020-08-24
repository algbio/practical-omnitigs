use crate::walks::VecWalk;
use traitgraph::index::GraphIndex;
use traitgraph::interface::{GraphBase, StaticGraph};

pub struct Unitigs<Graph: GraphBase> {
    unitigs: Vec<VecWalk<Graph>>,
}

impl<Graph: StaticGraph> Unitigs<Graph> {
    /// Computes the unitigs of a graph.
    /// Assumes that the graph is not a cycle, otherwise the method enters an endless loop.
    pub fn new(graph: &Graph) -> Self {
        let mut used_edges = vec![false; graph.edge_count()];
        let mut unitigs = Vec::new();

        for edge in graph.edge_indices() {
            if !used_edges[edge.as_usize()] {
                used_edges[edge.as_usize()] = true;

                let mut start_node = graph.edge_endpoints(edge).from_node;
                let mut end_node = graph.edge_endpoints(edge).to_node;
                let mut unitig = vec![end_node, start_node];

                while graph.is_biunivocal_node(start_node) && start_node != end_node {
                    let in_neighbor = graph.in_neighbors(start_node).into_iter().next().unwrap();
                    start_node = in_neighbor.node_id;
                    used_edges[in_neighbor.edge_id.as_usize()] = true;
                    unitig.push(start_node);
                }

                unitig.reverse();

                while graph.is_biunivocal_node(end_node) && start_node != end_node {
                    let out_neighbor = graph.out_neighbors(end_node).into_iter().next().unwrap();
                    end_node = out_neighbor.node_id;
                    used_edges[out_neighbor.edge_id.as_usize()] = true;
                    unitig.push(end_node);
                }

                unitigs.push(unitig.into());
            }
        }

        Self {unitigs}
    }
}

impl<'a, Graph: GraphBase> IntoIterator for &'a Unitigs<Graph> {
    type Item = &'a VecWalk<Graph>;
    type IntoIter = std::slice::Iter<'a, VecWalk<Graph>>;

    fn into_iter(self) -> Self::IntoIter {
        self.unitigs.iter()
    }
}
