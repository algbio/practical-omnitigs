use crate::interface::static_bigraph::StaticBigraph;
use crate::interface::static_bigraph::StaticEdgeCentricBigraph;
use crate::interface::static_bigraph::StaticNodeCentricBigraph;
use crate::interface::BidirectedData;
use traitgraph::interface::DynamicGraph;

pub trait DynamicBigraph: DynamicGraph + StaticBigraph {
    /// Make the nodes with the given two node ids mirror nodes.
    /// This may leave the old mirrors from a and b with dangling mirror pointers.
    fn set_mirror_nodes(&mut self, a: Self::NodeIndex, b: Self::NodeIndex);
}

pub trait DynamicNodeCentricBigraph: DynamicBigraph + StaticNodeCentricBigraph
where
    Self::NodeData: BidirectedData,
    Self::EdgeData: Clone,
{
    /**
     * Adds nodes such that the graph becomes a valid node-centric bigraph.
     * The indices of existing nodes are not altered.
     */
    fn add_mirror_nodes(&mut self) {
        for node_id in self.node_indices() {
            if self.mirror_node(node_id).is_none() {
                let mirror_index = self.add_node(self.node_data(node_id).reverse_complement());
                self.set_mirror_nodes(node_id, mirror_index);
            }
        }
    }

    /**
     * Adds edges such that the graph fulfils the node centric [mirror property].
     * Mirror edges get the cloned `EdgeData` from their existing mirror.
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn add_node_centric_mirror_edges(&mut self) {
        let mut edges = Vec::new();
        for from_id in self.node_indices() {
            let mut out_neighbors: Vec<_> = self.out_neighbors(from_id).into_iter().collect();
            out_neighbors.sort_by(|a, b| a.node_id.cmp(&b.node_id));
            out_neighbors.dedup_by(|a, b| a.node_id == b.node_id);
            for neighbor in out_neighbors {
                let to_id = neighbor.node_id;
                let mirror_from_id = self.mirror_node(to_id).unwrap();
                let mirror_to_id = self.mirror_node(from_id).unwrap();
                let difference = self
                    .edge_count_between(from_id, to_id)
                    .saturating_sub(self.edge_count_between(mirror_from_id, mirror_to_id));
                for _ in 0..difference {
                    edges.push((
                        mirror_from_id,
                        self.edge_data(neighbor.edge_id).clone(),
                        mirror_to_id,
                    ));
                }
            }
        }

        for edge in edges {
            self.add_edge(edge.0, edge.2, edge.1);
        }
    }
}

pub trait DynamicEdgeCentricBigraph: DynamicBigraph + StaticEdgeCentricBigraph
where
    Self::EdgeData: BidirectedData + Eq,
{
    /**
     * Adds edges such that the graph fulfils the edge centric [mirror property].
     * Mirror edges get the reverse complement `EdgeData` from their existing mirror.
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn add_edge_centric_mirror_edges(&mut self) {
        let mut edges = Vec::new();
        for from_id in self.node_indices() {
            for neighbor in self.out_neighbors(from_id) {
                let to_id = neighbor.node_id;
                let edge = neighbor.edge_id;
                if self.mirror_edge_edge_centric(edge).is_none() {
                    let mirror_from_id = self.mirror_node(to_id).unwrap();
                    let mirror_to_id = self.mirror_node(from_id).unwrap();
                    edges.push((
                        mirror_from_id,
                        self.edge_data(edge).reverse_complement(),
                        mirror_to_id,
                    ));
                }
            }
        }

        for edge in edges {
            self.add_edge(edge.0, edge.2, edge.1);
        }
    }
}

/*pub trait DynamicBigraphFromDigraph: DynamicBigraph + StaticBigraphFromDigraph + Sized
where
    Self::Topology: DynamicGraph,
    Self::EdgeData: Clone,
{
    /**
     * Converts the given topology into a bigraph.
     * Missing mirrors are added.
     */
    fn new_with_completed_nodes(
        _topology: Self::Topology,
        _binode_mapping_function: fn(&Self::NodeData) -> Self::NodeData,
    ) -> Self {
        unimplemented!()
    }

    /**
     * Converts the given topology into a bigraph that fulfils the [mirror property].
     * Missing mirrors and mirror edges are added.
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn new_with_completed_nodes_and_mirrored_edges(
        topology: Self::Topology,
        binode_mapping_function: fn(&Self::NodeData) -> Self::NodeData,
    ) -> Self {
        let mut bigraph = Self::new_with_completed_nodes(topology, binode_mapping_function);
        bigraph.add_node_centric_mirror_edges();
        bigraph
    }
}*/

#[cfg(test)]
mod tests {
    use crate::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
    use crate::interface::dynamic_bigraph::DynamicEdgeCentricBigraph;
    use crate::interface::dynamic_bigraph::DynamicNodeCentricBigraph;
    use crate::interface::static_bigraph::StaticEdgeCentricBigraph;
    use crate::interface::static_bigraph::StaticNodeCentricBigraph;
    use crate::interface::{
        static_bigraph::StaticBigraph, static_bigraph::StaticBigraphFromDigraph, BidirectedData,
    };
    use crate::traitgraph::implementation::petgraph_impl;
    use crate::traitgraph::interface::{ImmutableGraphContainer, MutableGraphContainer};

    #[test]
    fn test_bigraph_add_node_centric_mirror_edges() {
        let mut graph = petgraph_impl::new();
        #[derive(Eq, PartialEq, Debug, Hash, Clone)]
        struct NodeData(u32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
                Self(1000 - self.0)
            }
        }
        let n0 = graph.add_node(NodeData(0));
        let n1 = graph.add_node(NodeData(1));
        let n2 = graph.add_node(NodeData(2));
        let n3 = graph.add_node(NodeData(3));
        let n4 = graph.add_node(NodeData(997));
        graph.add_edge(n3, n4, ()); // This edge is a self-mirror
        graph.add_edge(n1, n2, ()); // This edge is not a self-mirror
        graph.add_edge(n1, n2, ()); // This edge is not a self-mirror
        graph.add_edge(n1, n2, ()); // This edge is not a self-mirror
        graph.add_edge(n0, n3, ()); // This edge is not a self-mirror
        graph.add_edge(n0, n3, ()); // This edge is not a self-mirror

        let mut graph = NodeBigraphWrapper::new_unchecked(graph);
        graph.add_mirror_nodes();
        assert!(graph.verify_node_pairing());
        assert_eq!(graph.node_count(), 8);

        graph.add_edge(
            graph.mirror_node(n2).unwrap(),
            graph.mirror_node(n1).unwrap(),
            (),
        );
        assert!(!graph.verify_node_mirror_property());
        graph.add_node_centric_mirror_edges();
        assert!(graph.verify_node_mirror_property());
        assert_eq!(graph.edge_count(), 11);
    }

    #[test]
    fn test_bigraph_add_edge_centric_mirror_edges() {
        let mut graph = petgraph_impl::new();
        #[derive(Eq, PartialEq, Debug, Hash, Clone)]
        struct NodeData(u32);
        impl BidirectedData for NodeData {
            fn reverse_complement(&self) -> Self {
                Self(1000 - self.0)
            }
        }
        #[derive(Eq, PartialEq, Debug, Hash, Clone)]
        struct EdgeData(u32);
        impl BidirectedData for EdgeData {
            fn reverse_complement(&self) -> Self {
                Self(1000 - self.0)
            }
        }
        let n0 = graph.add_node(NodeData(0));
        let n1 = graph.add_node(NodeData(1));
        let n2 = graph.add_node(NodeData(2));
        let n3 = graph.add_node(NodeData(3));
        let n4 = graph.add_node(NodeData(997));
        graph.add_edge(n3, n4, EdgeData(500)); // This edge is a self-mirror
        graph.add_edge(n1, n2, EdgeData(45)); // This edge is not a self-mirror
        graph.add_edge(n1, n2, EdgeData(194)); // This edge is not a self-mirror
        graph.add_edge(n1, n2, EdgeData(55)); // This edge is not a self-mirror
        graph.add_edge(n0, n3, EdgeData(2)); // This edge is not a self-mirror
        graph.add_edge(n0, n3, EdgeData(3)); // This edge is not a self-mirror

        let mut graph = NodeBigraphWrapper::new_unchecked(graph);
        graph.add_mirror_nodes();
        assert!(graph.verify_node_pairing());
        assert_eq!(graph.node_count(), 8);

        graph.add_edge(
            graph.mirror_node(n2).unwrap(),
            graph.mirror_node(n1).unwrap(),
            EdgeData(955),
        );
        assert!(!graph.verify_edge_mirror_property());
        graph.add_edge_centric_mirror_edges();
        assert!(graph.verify_edge_mirror_property());
        assert_eq!(graph.edge_count(), 12);
    }
}
