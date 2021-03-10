use crate::error::Result;
use bigraph::interface::dynamic_bigraph::{DynamicBigraph, DynamicEdgeCentricBigraph};
use bigraph::interface::BidirectedData;
use bigraph::traitgraph::algo::dijkstra::WeightedEdgeData;
use bigraph::traitgraph::interface::GraphBase;
use compact_genome::implementation::vector_genome_impl::VectorGenome;
use compact_genome::interface::Genome;
use std::collections::HashMap;
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader};
use std::iter::FromIterator;
use std::path::Path;
use std::rc::Rc;
//use bigraph::traitgraph::index::GraphIndex;

/// Type of graphs read from gfa files.
pub type PetGFAGraph<NodeData, EdgeData> =
    crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            BidirectedGFANodeData<NodeData>,
            EdgeData,
            usize,
        >,
    >;

/// The edge-centric variant of the type of graphs read from gfa files.
pub type PetGFAEdgeGraph<NodeData, EdgeData> =
    crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            NodeData,
            BidirectedGFANodeData<EdgeData>,
            usize,
        >,
    >;

/// Node data of a bidirected graph read from GFA
#[derive(Eq, PartialEq, Debug, Clone)]
pub struct BidirectedGFANodeData<T> {
    /// The sequence of this node. If forward is false, then this must be reverse complemented.
    pub sequence: Rc<VectorGenome>,
    /// True if this node is the forward node of sequence, false if it is the reverse complement node.
    pub forward: bool,
    /// Further data.
    pub data: T,
}

impl<T: BidirectedData> BidirectedData for BidirectedGFANodeData<T> {
    fn mirror(&self) -> Self {
        Self {
            sequence: self.sequence.clone(),
            forward: !self.forward,
            data: self.data.mirror(),
        }
    }
}

impl<T: WeightedEdgeData> WeightedEdgeData for BidirectedGFANodeData<T> {
    fn weight(&self) -> usize {
        self.data.weight()
    }
}

/// Read a bigraph in gfa format from a file.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_bigraph_from_file<
    P: AsRef<Path>,
    NodeData: Default,
    EdgeData: Default,
    Graph: DynamicBigraph<NodeData = BidirectedGFANodeData<NodeData>, EdgeData = EdgeData> + Default,
>(
    gfa_file: P,
) -> Result<(Graph, usize)> {
    read_gfa_as_bigraph(BufReader::new(File::open(gfa_file)?))
}

/// Read a bigraph in gfa format from a `BufRead`.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_bigraph<
    R: BufRead,
    NodeData: Default,
    EdgeData: Default,
    Graph: DynamicBigraph<NodeData = BidirectedGFANodeData<NodeData>, EdgeData = EdgeData> + Default,
>(
    gfa: R,
) -> Result<(Graph, usize)> {
    let mut graph = Graph::default();
    let mut k = usize::max_value();

    for line in gfa.lines() {
        let line = line?;

        if line.starts_with('H') {
            assert!(graph.is_empty());
            for column in line.split('\t') {
                if let Some(stripped) = column.strip_prefix("KL:Z:") {
                    assert_eq!(k, usize::max_value());
                    k = stripped.parse().unwrap();
                }
            }
        } else if line.starts_with('S') {
            assert_eq!(graph.edge_count(), 0);
            assert_ne!(k, usize::max_value());

            let mut columns = line.split('\t').skip(1);
            let node_index: usize = columns.next().unwrap().parse().unwrap();
            assert_eq!((node_index - 1) * 2, graph.node_count());

            let sequence = columns.next().unwrap();
            let sequence = Rc::new(VectorGenome::from_iter(sequence.bytes()));
            assert!(columns.next().is_none());
            assert!(
                sequence.len() >= k,
                "Node {} has sequence '{}' of length {} (k = {})",
                node_index,
                sequence,
                sequence.len(),
                k
            );

            let n1 = graph.add_node(BidirectedGFANodeData {
                sequence: sequence.clone(),
                forward: true,
                data: Default::default(),
            });
            let n2 = graph.add_node(BidirectedGFANodeData {
                sequence: sequence.clone(),
                forward: false,
                data: Default::default(),
            });
            graph.set_mirror_nodes(n1, n2);
        } else if line.starts_with('L') {
            assert_ne!(k, usize::max_value());

            let mut columns = line.split('\t').skip(1);
            let n1 = (columns.next().unwrap().parse::<usize>().unwrap() * 2
                - if columns.next().unwrap() == "+" { 2 } else { 1 })
            .into();
            let n2 = (columns.next().unwrap().parse::<usize>().unwrap() * 2
                - if columns.next().unwrap() == "+" { 2 } else { 1 })
            .into();

            let has_edge = graph.contains_edge_between(n1, n2);
            assert_eq!(
                has_edge,
                graph.contains_edge_between(
                    graph.mirror_node(n2).unwrap(),
                    graph.mirror_node(n1).unwrap()
                )
            );

            if !has_edge {
                graph.add_edge(n1, n2, Default::default());
                graph.add_edge(
                    graph.mirror_node(n2).unwrap(),
                    graph.mirror_node(n1).unwrap(),
                    Default::default(),
                );
            }
        }
    }

    Ok((graph, k))
}

/// Read an edge-centric bigraph in gfa format from a file.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_edge_centric_bigraph_from_file<
    P: AsRef<Path>,
    NodeData: Default,
    EdgeData: Default + BidirectedData + Eq + Clone,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = BidirectedGFANodeData<EdgeData>>
        + Default
        + std::fmt::Debug,
>(
    gfa_file: P,
) -> Result<(Graph, usize)> {
    read_gfa_as_edge_centric_bigraph(BufReader::new(File::open(gfa_file)?))
}

fn get_or_create_node<Graph: DynamicBigraph, G: Genome + Hash>(
    bigraph: &mut Graph,
    id_map: &mut HashMap<G, <Graph as GraphBase>::NodeIndex>,
    genome: &G,
) -> <Graph as GraphBase>::NodeIndex
where
    for<'a> &'a G: IntoIterator<Item = u8>,
    <Graph as GraphBase>::NodeData: Default,
    <Graph as GraphBase>::EdgeData: Clone,
{
    if let Some(node) = id_map.get(genome) {
        *node
    } else {
        let node = bigraph.add_node(Default::default());
        id_map.insert(genome.clone(), node);

        let reverse_complement = genome.reverse_complement();
        if &reverse_complement == genome {
            bigraph.set_mirror_nodes(node, node);
        } else {
            let mirror_node = bigraph.add_node(Default::default());
            id_map.insert(reverse_complement, mirror_node);
            bigraph.set_mirror_nodes(node, mirror_node);
        }

        node
    }
}

/// Read an edge-centric bigraph in gfa format from a `BufRead`.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_edge_centric_bigraph<
    R: BufRead,
    NodeData: Default,
    EdgeData: Default + BidirectedData + Eq + Clone,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = BidirectedGFANodeData<EdgeData>>
        + Default
        + std::fmt::Debug,
>(
    gfa: R,
) -> Result<(Graph, usize)> {
    let mut bigraph = Graph::default();
    let mut id_map = HashMap::new();
    let mut k = usize::max_value();

    for line in gfa.lines() {
        let line = line?;

        if line.starts_with('H') {
            assert!(bigraph.is_empty());
            for column in line.split('\t') {
                if let Some(stripped) = column.strip_prefix("KL:Z:") {
                    assert_eq!(k, usize::max_value());
                    k = stripped.parse().unwrap();
                }
            }
        } else if line.starts_with('S') {
            assert_ne!(k, usize::max_value());

            let mut columns = line.split('\t').skip(1);
            let node_index: usize = columns.next().unwrap().parse().unwrap();
            assert_eq!((node_index - 1) * 2, bigraph.edge_count());

            let sequence = columns.next().unwrap();
            //println!("sequence {}", sequence);
            let sequence = Rc::new(VectorGenome::from_iter(sequence.bytes()));
            let edge_data = BidirectedGFANodeData {
                sequence: sequence.clone(),
                forward: true,
                data: Default::default(),
            };
            let reverse_edge_data = edge_data.mirror();

            assert!(columns.next().is_none());
            assert!(
                sequence.len() >= k,
                "Node {} has sequence '{}' of length {} (k = {})",
                node_index,
                sequence,
                sequence.len(),
                k
            );

            let pre_plus = sequence.prefix(k - 1);
            let pre_minus = sequence.suffix(k - 1).reverse_complement();
            let succ_plus = sequence.suffix(k - 1);
            let succ_minus = sequence.prefix(k - 1).reverse_complement();

            let pre_plus = get_or_create_node(&mut bigraph, &mut id_map, &pre_plus);
            let pre_minus = get_or_create_node(&mut bigraph, &mut id_map, &pre_minus);
            let succ_plus = get_or_create_node(&mut bigraph, &mut id_map, &succ_plus);
            let succ_minus = get_or_create_node(&mut bigraph, &mut id_map, &succ_minus);

            //println!("Adding edge ({}, {}) and reverse ({}, {})", pre_plus.as_usize(), succ_plus.as_usize(), pre_minus.as_usize(), succ_minus.as_usize());
            bigraph.add_edge(pre_plus, succ_plus, edge_data);
            bigraph.add_edge(pre_minus, succ_minus, reverse_edge_data);
        } else if line.starts_with('L') {
            assert_ne!(k, usize::max_value());

            // Since we are using a hashtable to find the nodes, we can ignore the edges.
        }
    }

    //println!("{:?}", bigraph);
    assert!(bigraph.verify_node_pairing());
    assert!(bigraph.verify_edge_mirror_property());
    Ok((bigraph, k))
}

#[cfg(test)]
mod tests {
    use crate::io::gfa::{read_gfa_as_edge_centric_bigraph, PetGFAEdgeGraph};
    use std::io::BufReader;

    #[test]
    fn test_read_gfa_as_edge_centric_bigraph_simple() {
        let gfa = "H\tKL:Z:3\nS\t1\tACGA\nS\t2\tTCGT";
        let (_bigraph, k): (PetGFAEdgeGraph<(), ()>, _) =
            read_gfa_as_edge_centric_bigraph(BufReader::new(gfa.as_bytes())).unwrap();
        assert_eq!(k, 3);
    }
}
