use traitgraph::index::GraphIndex;
use traitgraph::interface::{GraphBase, StaticGraph};
use traitgraph::walks::{NodeWalk, VecNodeWalk};

/// Algorithms to count uncompacted unitigs in a graph.
pub mod uncompacted_unitigs;

/// A unitig stored as sequence of nodes.
#[derive(Debug, Clone, Sequence)]
pub struct NodeUnitig<Graph: GraphBase> {
    /// Store the unitig as a node-centric walk.
    /// This allows unitigs that are a single node to be stored.
    #[inner_sequence]
    walk: VecNodeWalk<Graph>,
    /// If a unitig has just a single edge, this edge might be a multiedge.
    /// In this case, this field is required to uniquely identify the unitig.
    single_edge_disambiguator: Option<Graph::EdgeIndex>,
}

impl<Graph: GraphBase> NodeUnitig<Graph> {
    /// Creates a new unitig from the given walk.
    /// Panics if the walk as length two nodes but no single edge is given, or if a single edge is given but the walk has a length other than two nodes.
    pub fn new(walk: VecNodeWalk<Graph>, single_edge: Option<Graph::EdgeIndex>) -> Self {
        assert_eq!(walk.len() == 2, single_edge.is_some());
        Self {
            walk,
            single_edge_disambiguator: single_edge,
        }
    }

    /// Returns an iterator over the nodes of this unitig.
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = Graph::NodeIndex> {
        self.walk.iter()
    }

    /// Returns the single edge field of this unitig.
    /// The field is `Some` if this unitig has exactly two nodes.
    /// In this case, the two nodes might be connected by multiple edges, so then this field serves as disambiguation between the edges.
    pub fn single_edge(&self) -> &Option<Graph::EdgeIndex> {
        &self.single_edge_disambiguator
    }

    /// Extracts the node walk from this unitig, consuming the unitig.
    /// Note that the information for unitigs of length two with a multiedge is lost.
    pub fn into_node_walk(self) -> VecNodeWalk<Graph> {
        self.walk
    }

    /// Returns the length of this unitigs as its amount of nodes.
    pub fn len(&self) -> usize {
        self.walk.len()
    }

    /// Returns true if this unitig contains no nodes.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A structure storing a vector of node-centric unitigs.
pub struct Unitigs<Graph: GraphBase> {
    unitigs: Vec<NodeUnitig<Graph>>,
}

impl<Graph: StaticGraph> Unitigs<Graph> {
    /// Computes the maximal unitigs of a graph.
    ///
    /// The unitigs are computed both node- and edge-centric.
    /// That means that parallel edges are counted as separate unitig each, and single nodes are counted as unitigs as well.
    ///
    /// If the graph is a cycle, it outputs an arbitrary subwalk.
    pub fn new(graph: &Graph) -> Self {
        let mut used_edges = vec![false; graph.edge_count()];
        let mut unitigs = Vec::new();

        // Add unitigs of length at least two nodes.
        for edge in graph.edge_indices() {
            if !used_edges[edge.as_usize()] {
                used_edges[edge.as_usize()] = true;

                let mut start_node = graph.edge_endpoints(edge).from_node;
                let mut end_node = graph.edge_endpoints(edge).to_node;
                let mut unitig = vec![end_node, start_node];

                while graph.is_biunivocal_node(start_node) && start_node != end_node {
                    let in_neighbor = graph.in_neighbors(start_node).next().unwrap();
                    start_node = in_neighbor.node_id;
                    used_edges[in_neighbor.edge_id.as_usize()] = true;
                    unitig.push(start_node);
                }

                unitig.reverse();

                while graph.is_biunivocal_node(end_node) && start_node != end_node {
                    let out_neighbor = graph.out_neighbors(end_node).next().unwrap();
                    end_node = out_neighbor.node_id;
                    used_edges[out_neighbor.edge_id.as_usize()] = true;
                    unitig.push(end_node);
                }

                let single_edge_disambiguator = if unitig.len() == 2 { Some(edge) } else { None };
                unitigs.push(NodeUnitig::new(unitig.into(), single_edge_disambiguator));
            }
        }

        // Add single nodes
        for node in graph.node_indices() {
            if graph.is_bivalent_node(node) {
                unitigs.push(NodeUnitig::new(vec![node].into(), None));
            }
        }

        Self { unitigs }
    }

    /// Returns an iterator over the nodes of this unitig.
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = &'a NodeUnitig<Graph>> {
        self.unitigs.iter()
    }
}

impl<Graph: GraphBase> IntoIterator for Unitigs<Graph> {
    type Item = NodeUnitig<Graph>;
    type IntoIter = std::vec::IntoIter<NodeUnitig<Graph>>;

    fn into_iter(self) -> Self::IntoIter {
        self.unitigs.into_iter()
    }
}

impl<Graph: GraphBase> PartialEq for NodeUnitig<Graph>
where
    Graph::NodeIndex: PartialEq,
    Graph::EdgeIndex: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.walk == rhs.walk && self.single_edge_disambiguator == rhs.single_edge_disambiguator
    }
}

impl<Graph: GraphBase> Eq for NodeUnitig<Graph>
where
    Graph::NodeIndex: Eq,
    Graph::EdgeIndex: Eq,
{
}

#[cfg(test)]
mod test {
    use super::{NodeUnitig, Unitigs};
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
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n3, n5, 14);
        graph.add_edge(n4, n8, 15);
        graph.add_edge(n5, n8, 16);
        let e7 = graph.add_edge(n8, n6, 17);
        let e8 = graph.add_edge(n8, n6, 175);
        graph.add_edge(n8, n7, 18);
        let e10 = graph.add_edge(n6, n0, 19);
        graph.add_edge(n7, n0, 20);

        let unitigs = Unitigs::new(&graph);
        let mut unitigs_iter = unitigs.iter();
        assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n0, n1, n2, n3]),
                None
            ))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n3, n4, n8]),
                None
            ))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n3, n5, n8]),
                None
            ))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n8, n6]),
                Some(e7)
            ))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n8, n6]),
                Some(e8)
            ))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n8, n7, n0]),
                None
            ))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(
                graph.create_node_walk(&[n6, n0]),
                Some(e10)
            ))
        );
        assert_eq!(
            unitigs_iter.next(),
            Some(&NodeUnitig::new(graph.create_node_walk(&[n8]), None))
        );
        assert_eq!(unitigs_iter.next(), None);
    }
}
