use super::{MacronodeAlgorithm, Macronodes};
use crate::unitigs::{NodeUnitig, NodeUnitigs};
use traitgraph::interface::StaticGraph;
use traitsequence::interface::Sequence;

/// Compute the macronodes of a strongly connected graph.
pub struct StronglyConnectedMacronodes;

impl<Graph: StaticGraph> MacronodeAlgorithm<Graph> for StronglyConnectedMacronodes {
    fn compute_macronodes(graph: &Graph) -> Macronodes<Graph> {
        let unitigs = NodeUnitigs::compute(graph);
        let macronodes: Vec<_> = unitigs
            .into_iter()
            .filter(|unitig| {
                (graph.out_degree(*unitig.iter().next().unwrap()) == 1
                    && graph.in_degree(*unitig.iter().last().unwrap()) == 1)
                    || (unitig.len() == 1 && graph.is_bivalent_node(*unitig.iter().next().unwrap()))
            })
            .map(NodeUnitig::into_node_walk)
            .collect();
        Macronodes::from(macronodes)
    }
}

#[cfg(test)]
mod tests {
    use super::StronglyConnectedMacronodes;
    use crate::macrotigs::macronodes::MacronodeAlgorithm;
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::{MutableGraphContainer, WalkableGraph};
    use traitsequence::interface::Sequence;

    #[test]
    fn test_compute_macronodes_complex_graph() {
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
        let n9 = graph.add_node(9);
        let n10 = graph.add_node(10);
        let n11 = graph.add_node(11);
        let n12 = graph.add_node(12);
        let n13 = graph.add_node(13);
        let n14 = graph.add_node(14);
        let n15 = graph.add_node(15);
        let n16 = graph.add_node(16);
        let n17 = graph.add_node(17);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n4, 13);
        graph.add_edge(n3, n5, 14);
        graph.add_edge(n4, n8, 15);
        graph.add_edge(n5, n8, 16);
        graph.add_edge(n8, n6, 17);
        graph.add_edge(n8, n6, 175);
        graph.add_edge(n8, n7, 18);
        graph.add_edge(n6, n0, 19);
        graph.add_edge(n7, n0, 20);
        graph.add_edge(n8, n9, 21);
        graph.add_edge(n9, n10, 22);
        graph.add_edge(n10, n8, 23);
        graph.add_edge(n11, n4, 24);
        graph.add_edge(n11, n5, 25);
        graph.add_edge(n6, n11, 26);
        graph.add_edge(n7, n11, 27);
        graph.add_edge(n8, n12, 28);
        graph.add_edge(n8, n12, 29);
        graph.add_edge(n12, n13, 30);
        graph.add_edge(n13, n14, 31);
        graph.add_edge(n14, n8, 32);
        graph.add_edge(n8, n15, 33);
        graph.add_edge(n15, n16, 34);
        graph.add_edge(n16, n17, 35);
        graph.add_edge(n17, n8, 36);
        graph.add_edge(n17, n8, 37);

        let macronodes = StronglyConnectedMacronodes::compute_macronodes(&graph);
        let mut macronodes_iter = macronodes.iter();
        debug_assert_eq!(
            macronodes_iter.next(),
            Some(&graph.create_node_walk(&[n0, n1, n2, n3]))
        );
        debug_assert_eq!(macronodes_iter.next(), Some(&graph.create_node_walk(&[n6])));
        debug_assert_eq!(macronodes_iter.next(), Some(&graph.create_node_walk(&[n8])));
        debug_assert_eq!(
            macronodes_iter.next(),
            Some(&graph.create_node_walk(&[n11]))
        );
        debug_assert_eq!(macronodes_iter.next(), None);
    }

    #[test]
    fn test_compute_macronodes_empty_graph() {
        let graph = petgraph_impl::new::<i32, i32>();

        let macronodes = StronglyConnectedMacronodes::compute_macronodes(&graph);
        let mut macronodes_iter = macronodes.iter();
        debug_assert_eq!(macronodes_iter.next(), None);
    }
}
