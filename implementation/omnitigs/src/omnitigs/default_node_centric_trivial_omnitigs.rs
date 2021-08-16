use crate::omnitigs::univocal_extension_algorithms::{
    NonSccNodeCentricUnivocalExtensionStrategy, SccNodeCentricUnivocalExtensionStrategy,
};
use crate::omnitigs::{NodeCentricUnivocalExtensionAlgorithm, TrivialNodeCentricOmnitigAlgorithm};
use bitvector::BitVector;
use std::marker::PhantomData;
use traitgraph::index::GraphIndex;
use traitgraph::interface::StaticGraph;
use traitgraph::walks::VecNodeWalk;
use traitsequence::interface::Sequence;

/// An algorithm to extract trivial omnitigs.
pub struct DefaultTrivialNodeCentricOmnitigAlgorithm<NodeCentricUnivocalExtensionStrategy> {
    _univocal_extension_strategy: PhantomData<NodeCentricUnivocalExtensionStrategy>,
}

/// An algorithm to extract trivial omnitigs form a strongly connected graph.
pub type SccTrivialOmnitigAlgorithm =
    DefaultTrivialNodeCentricOmnitigAlgorithm<SccNodeCentricUnivocalExtensionStrategy>;

/// An algorithm to extract trivial omnitigs form a not strongly connected graph.
/// This runs slightly slower than the counterpart for strongly connected graphs, especially for long univocal extensions.
pub type NonSccTrivialOmnitigAlgorithm =
    DefaultTrivialNodeCentricOmnitigAlgorithm<NonSccNodeCentricUnivocalExtensionStrategy>;

impl<
        Graph: StaticGraph,
        NodeCentricUnivocalExtensionStrategy: NodeCentricUnivocalExtensionAlgorithm<Graph, VecNodeWalk<Graph>>,
    > TrivialNodeCentricOmnitigAlgorithm<Graph>
    for DefaultTrivialNodeCentricOmnitigAlgorithm<NodeCentricUnivocalExtensionStrategy>
{
    type NodeCentricUnivocalExtensionStrategy = NodeCentricUnivocalExtensionStrategy;

    fn compute_maximal_trivial_node_centric_omnitigs(
        graph: &Graph,
        mut omnitigs: Vec<VecNodeWalk<Graph>>,
    ) -> Vec<VecNodeWalk<Graph>> {
        info!("Marking used edges");
        let mut used_nodes = BitVector::new(graph.node_count());
        for omnitig in omnitigs.iter() {
            for node in omnitig.iter() {
                used_nodes.insert(node.as_usize());
            }
        }

        info!("Extend {} unused nodes", graph.node_count());
        for node in graph.node_indices() {
            if used_nodes.contains(node.as_usize()) {
                continue;
            }

            let mut trivial_omnitig: VecNodeWalk<Graph> =
                Self::NodeCentricUnivocalExtensionStrategy::compute_univocal_extension(
                    graph,
                    &[node],
                );

            let is_univocal = trivial_omnitig
                .iter()
                .rev()
                .skip(1)
                .all(|&n| graph.out_degree(n) == 1);

            if is_univocal {
                let prefix = &mut trivial_omnitig;
                prefix.reverse();

                let mut can_extend = true;
                while can_extend {
                    can_extend = false;
                    for in_neighbor in graph.in_neighbors(*prefix.last().unwrap()) {
                        if graph.out_degree(in_neighbor.node_id) == 1
                            && !prefix.contains(&in_neighbor.node_id)
                        {
                            prefix.push(in_neighbor.node_id);
                            can_extend = true;
                            break;
                        }
                    }
                }

                prefix.reverse();
            }

            let is_r_univocal = trivial_omnitig
                .iter()
                .skip(1)
                .all(|&n| graph.in_degree(n) == 1);

            if is_r_univocal {
                let suffix = &mut trivial_omnitig;

                let mut can_extend = true;
                while can_extend {
                    can_extend = false;
                    for out_neighbor in graph.out_neighbors(*suffix.last().unwrap()) {
                        if graph.in_degree(out_neighbor.node_id) == 1
                            && !suffix.contains(&out_neighbor.node_id)
                        {
                            suffix.push(out_neighbor.node_id);
                            can_extend = true;
                            break;
                        }
                    }
                }
            }

            /*let last_split_edge = trivial_omnitig
                .iter()
                .enumerate()
                .filter(|(_, e)| graph.is_split_edge(**e))
                .map(|(i, _)| i)
                .last()
                .unwrap_or(0);
            let first_join_edge = trivial_omnitig
                .iter()
                .enumerate()
                .filter(|(_, e)| graph.is_join_edge(**e))
                .map(|(i, _)| i)
                .next()
                .unwrap_or(trivial_omnitig.len() - 1);*/

            for node in trivial_omnitig.iter() {
                used_nodes.insert(node.as_usize());
            }
            omnitigs.push(trivial_omnitig);
        }

        debug_assert_eq!(used_nodes.len(), graph.node_count());

        omnitigs
    }
}

#[cfg(test)]
mod tests {
    use crate::omnitigs::NodeCentricOmnitigs;
    use traitgraph::implementation::petgraph_impl;
    use traitgraph::interface::MutableGraphContainer;

    #[test]
    fn test_compute_trivial_omnitigs_cycle() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n1, n2, ());
        graph.add_edge(n2, n0, ());

        let trivial_omnitigs = Vec::compute_trivial_node_centric_omnitigs(&graph);
        debug_assert_eq!(trivial_omnitigs, vec![vec![n1, n2, n0, n1, n2]]);
    }

    #[test]
    fn test_compute_trivial_omnitigs_trap_cycle() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        graph.add_edge(n3, n0, ());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n1, n2, ());
        graph.add_edge(n2, n0, ());

        let trivial_omnitigs = Vec::compute_trivial_node_centric_omnitigs_non_scc(&graph);
        debug_assert_eq!(trivial_omnitigs, vec![vec![n3, n0, n1, n2]]);
    }

    #[test]
    fn test_compute_trivial_omnitigs_reverse_trap_cycle() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        graph.add_edge(n0, n3, ());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n1, n2, ());
        graph.add_edge(n2, n0, ());

        let trivial_omnitigs = Vec::compute_trivial_node_centric_omnitigs_non_scc(&graph);
        debug_assert_eq!(trivial_omnitigs, vec![vec![n1, n2, n0, n3]]);
    }

    #[test]
    fn test_compute_trivial_omnitigs_trap_cycle_inverse_id_order() {
        let mut graph = petgraph_impl::new();
        let n3 = graph.add_node(());
        let n2 = graph.add_node(());
        let n1 = graph.add_node(());
        let n0 = graph.add_node(());
        graph.add_edge(n2, n0, ());
        graph.add_edge(n1, n2, ());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n3, n0, ());

        let trivial_omnitigs = Vec::compute_trivial_node_centric_omnitigs_non_scc(&graph);
        debug_assert_eq!(trivial_omnitigs, vec![vec![n3, n0, n1, n2]]);
    }

    #[test]
    fn test_compute_trivial_omnitigs_reverse_trap_cycle_inverse_id_order() {
        let mut graph = petgraph_impl::new();
        let n3 = graph.add_node(());
        let n2 = graph.add_node(());
        let n1 = graph.add_node(());
        let n0 = graph.add_node(());
        graph.add_edge(n2, n0, ());
        graph.add_edge(n1, n2, ());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n0, n3, ());

        let trivial_omnitigs = Vec::compute_trivial_node_centric_omnitigs_non_scc(&graph);
        debug_assert_eq!(trivial_omnitigs, vec![vec![n1, n2, n0, n3]]);
    }

    #[test]
    fn test_compute_trivial_omnitigs_path() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n1, n2, ());

        let trivial_omnitigs = Vec::compute_trivial_node_centric_omnitigs_non_scc(&graph);
        debug_assert_eq!(trivial_omnitigs, vec![vec![n0, n1, n2]]);
    }
}
