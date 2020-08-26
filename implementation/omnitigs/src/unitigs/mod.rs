use traitgraph::index::GraphIndex;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::VecEdgeWalk;

pub mod uncompacted_unitigs;

pub type Unitig<Graph> = VecEdgeWalk<Graph>;

pub struct Unitigs<Graph: GraphBase> {
    unitigs: Vec<VecEdgeWalk<Graph>>,
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
                let mut unitig = vec![edge];

                while graph.is_biunivocal_node(start_node) && start_node != end_node {
                    let in_neighbor = graph.in_neighbors(start_node).into_iter().next().unwrap();
                    start_node = in_neighbor.node_id;
                    used_edges[in_neighbor.edge_id.as_usize()] = true;
                    unitig.push(in_neighbor.edge_id);
                }

                unitig.reverse();

                while graph.is_biunivocal_node(end_node) && start_node != end_node {
                    let out_neighbor = graph.out_neighbors(end_node).into_iter().next().unwrap();
                    end_node = out_neighbor.node_id;
                    used_edges[out_neighbor.edge_id.as_usize()] = true;
                    unitig.push(out_neighbor.edge_id);
                }

                unitigs.push(unitig.into());
            }
        }

        Self { unitigs }
    }

    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = &'a VecEdgeWalk<Graph>> {
        self.unitigs.iter()
    }
}

impl<Graph: GraphBase> IntoIterator for Unitigs<Graph> {
    type Item = VecEdgeWalk<Graph>;
    type IntoIter = std::vec::IntoIter<VecEdgeWalk<Graph>>;

    fn into_iter(self) -> Self::IntoIter {
        self.unitigs.into_iter()
    }
}

#[cfg(test)]
mod test {
    use super::Unitigs;
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::{MutableGraphContainer, WalkableGraph};

    #[test]
    fn test_unitig_computation() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        let n4 = graph.add_node(4);
        let n5 = graph.add_node(5);
        let n6 = graph.add_node(6);
        let n7 = graph.add_node(7);
        let n8 = graph.add_node(8);
        let e0 = graph.add_edge(n0, n1, 10);
        let e1 = graph.add_edge(n1, n2, 11);
        let e2 = graph.add_edge(n2, n3, 12);
        let e3 = graph.add_edge(n3, n4, 13);
        let e4 = graph.add_edge(n3, n5, 14);
        let e5 = graph.add_edge(n4, n8, 15);
        let e6 = graph.add_edge(n5, n8, 16);
        let e7 = graph.add_edge(n8, n6, 17);
        let e8 = graph.add_edge(n8, n6, 175);
        let e9 = graph.add_edge(n8, n7, 18);
        let e10 = graph.add_edge(n6, n0, 19);
        let e11 = graph.add_edge(n7, n0, 20);

        let unitigs = Unitigs::new(&graph);
        let mut unitigs_iter = unitigs.iter();
        // TODO move walks to traitgraph and support properly from another graph trait, i.e. WalkableGraph.
        assert_eq!(
            unitigs_iter.next(),
            Some(&graph.create_edge_walk(&[e0, e1, e2]))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&graph.create_edge_walk(&[e3, e5]))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&graph.create_edge_walk(&[e4, e6]))
        );
        assert_eq!(unitigs_iter.next(), Some(&graph.create_edge_walk(&[e7])));
        assert_eq!(unitigs_iter.next(), Some(&graph.create_edge_walk(&[e8])));
        assert_eq!(
            unitigs_iter.next(),
            Some(&graph.create_edge_walk(&[e9, e11]))
        );
        assert_eq!(unitigs_iter.next(), Some(&graph.create_edge_walk(&[e10])));
        assert_eq!(unitigs_iter.next(), None);
    }
}
