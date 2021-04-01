use crate::error::Result;
use crate::io::fasta::FastaData;
use bigraph::interface::dynamic_bigraph::{DynamicBigraph, DynamicEdgeCentricBigraph};
use bigraph::interface::BidirectedData;
use bigraph::traitgraph::algo::dijkstra::WeightedEdgeData;
use bigraph::traitgraph::index::GraphIndex;
use bigraph::traitgraph::interface::GraphBase;
use bigraph::traitgraph::traitsequence::interface::Sequence;
use compact_genome::implementation::vec_sequence::AsciiVectorGenome;
use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
use std::collections::HashMap;
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::rc::Rc;
//use bigraph::traitgraph::index::GraphIndex;

/// Type of graphs read from gfa files.
pub type PetGfaGraph<NodeData, EdgeData> =
    crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            BidirectedGfaNodeData<NodeData>,
            EdgeData,
            usize,
        >,
    >;

/// The edge-centric variant of the type of graphs read from gfa files.
pub type PetGfaEdgeGraph<NodeData, EdgeData> =
    crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            NodeData,
            BidirectedGfaNodeData<EdgeData>,
            usize,
        >,
    >;

/// Node data of a bidirected graph read from GFA
#[derive(Eq, PartialEq, Debug, Clone, Default)]
pub struct BidirectedGfaNodeData<T> {
    /// The sequence of this node. If forward is false, then this must be reverse complemented.
    pub sequence: Rc<AsciiVectorGenome>,
    /// True if this node is the forward node of sequence, false if it is the reverse complement node.
    pub forward: bool,
    /// Further data.
    pub data: T,
}

impl<T: BidirectedData> BidirectedData for BidirectedGfaNodeData<T> {
    fn mirror(&self) -> Self {
        Self {
            sequence: self.sequence.clone(),
            forward: !self.forward,
            data: self.data.mirror(),
        }
    }
}

impl<T: WeightedEdgeData> WeightedEdgeData for BidirectedGfaNodeData<T> {
    fn weight(&self) -> usize {
        self.data.weight()
    }
}

impl<T> FastaData for BidirectedGfaNodeData<T> {
    type Genome = AsciiVectorGenome;
    type GenomeSubsequence = [u8];

    fn sequence(&self) -> &Self::Genome {
        &self.sequence
    }
}

/// Read a bigraph in gfa format from a file.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_bigraph_from_file<
    P: AsRef<Path>,
    NodeData: Default,
    EdgeData: Default,
    Graph: DynamicBigraph<NodeData = BidirectedGfaNodeData<NodeData>, EdgeData = EdgeData> + Default,
>(
    gfa_file: P,
    ignore_k: bool,
    allow_messy_edges: bool,
) -> Result<(Graph, usize)> {
    read_gfa_as_bigraph(
        BufReader::new(File::open(gfa_file)?),
        ignore_k,
        allow_messy_edges,
    )
}

/// Read a bigraph in gfa format from a `BufRead`.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_bigraph<
    R: BufRead,
    NodeData: Default,
    EdgeData: Default,
    Graph: DynamicBigraph<NodeData = BidirectedGfaNodeData<NodeData>, EdgeData = EdgeData> + Default,
>(
    gfa: R,
    ignore_k: bool,
    allow_messy_edges: bool,
) -> Result<(Graph, usize)> {
    let mut graph = Graph::default();
    let mut k = usize::max_value();
    let mut node_name_map = HashMap::new();

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
            if !allow_messy_edges {
                assert_eq!(graph.edge_count(), 0);
            }
            if !ignore_k {
                assert_ne!(k, usize::max_value());
            }

            let mut columns = line.split('\t').skip(1);
            let node_name: &str = columns.next().unwrap();

            let sequence = columns.next().unwrap();
            let sequence: AsciiVectorGenome = sequence.bytes().collect();
            let sequence = Rc::new(sequence);
            assert!(
                sequence.len() >= k || ignore_k,
                "Node {} has sequence '{:?}' of length {} (k = {})",
                node_name,
                sequence,
                sequence.len(),
                k
            );

            let n1 = graph.add_node(BidirectedGfaNodeData {
                sequence: sequence.clone(),
                forward: true,
                data: Default::default(),
            });
            let n2 = graph.add_node(BidirectedGfaNodeData {
                sequence: sequence.clone(),
                forward: false,
                data: Default::default(),
            });
            graph.set_mirror_nodes(n1, n2);
            node_name_map.insert(node_name.to_owned(), n1);
        } else if line.starts_with('L') {
            if !ignore_k {
                assert_ne!(k, usize::max_value());
            }

            let mut columns = line.split('\t').skip(1);
            let n1_name = columns.next().unwrap();
            let n1_direction = if columns.next().unwrap() == "+" { 0 } else { 1 };
            let n2_name = columns.next().unwrap();
            let n2_direction = if columns.next().unwrap() == "+" { 0 } else { 1 };

            if let (Some(n1), Some(n2)) = (node_name_map.get(n1_name), node_name_map.get(n2_name)) {
                let n1 = (n1.as_usize() + n1_direction).into();
                let n2 = (n2.as_usize() + n2_direction).into();

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
    }

    if ignore_k {
        k = 0;
    }

    Ok((graph, k))
}

/// Read an edge-centric bigraph in gfa format from a file.
/// This method also returns the k-mer length given in the gfa file as well as the full gfa header.
pub fn read_gfa_as_edge_centric_bigraph_from_file<
    P: AsRef<Path>,
    NodeData: Default,
    EdgeData: Default + BidirectedData + Eq + Clone,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = BidirectedGfaNodeData<EdgeData>>
        + Default
        + std::fmt::Debug,
>(
    gfa_file: P,
    estimate_k: bool,
) -> Result<(Graph, usize, String)> {
    read_gfa_as_edge_centric_bigraph(BufReader::new(File::open(gfa_file)?), estimate_k)
}

fn get_or_create_node<
    Graph: DynamicBigraph,
    G: for<'a> OwnedGenomeSequence<'a, GenomeSubsequence> + Hash + Eq + Clone,
    GenomeSubsequence: for<'a> GenomeSequence<'a, GenomeSubsequence> + ?Sized,
>(
    bigraph: &mut Graph,
    id_map: &mut HashMap<G, <Graph as GraphBase>::NodeIndex>,
    genome: G,
) -> <Graph as GraphBase>::NodeIndex
where
    <Graph as GraphBase>::NodeData: Default,
    <Graph as GraphBase>::EdgeData: Clone,
{
    if let Some(node) = id_map.get(&genome) {
        *node
    } else {
        let node = bigraph.add_node(Default::default());

        let reverse_complement = genome.reverse_complement();
        if reverse_complement == genome {
            bigraph.set_mirror_nodes(node, node);
        } else {
            let mirror_node = bigraph.add_node(Default::default());
            id_map.insert(reverse_complement, mirror_node);
            bigraph.set_mirror_nodes(node, mirror_node);
        }

        id_map.insert(genome, node);

        node
    }
}

/// Read an edge-centric bigraph in gfa format from a `BufRead`.
/// This method also returns the k-mer length given in the gfa file as well as the full gfa header.
pub fn read_gfa_as_edge_centric_bigraph<
    R: BufRead,
    NodeData: Default,
    EdgeData: Default + BidirectedData + Eq + Clone,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = BidirectedGfaNodeData<EdgeData>>
        + Default
        + std::fmt::Debug,
>(
    gfa: R,
    estimate_k: bool,
) -> Result<(Graph, usize, String)> {
    assert!(!estimate_k, "Estimating k not supported yet");

    let mut bigraph = Graph::default();
    let mut id_map = HashMap::new();
    let mut k = usize::max_value();
    let mut header = None;

    for line in gfa.lines() {
        let line = line?;

        if line.starts_with('H') {
            assert!(bigraph.is_empty());
            header = Some(line.clone());
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
            let sequence: AsciiVectorGenome = sequence.bytes().collect();
            let sequence = Rc::new(sequence);
            let edge_data = BidirectedGfaNodeData {
                sequence: sequence.clone(),
                forward: true,
                data: Default::default(),
            };
            let reverse_edge_data = edge_data.mirror();

            assert!(columns.next().is_none());
            assert!(
                sequence.len() >= k,
                "Node {} has sequence '{:?}' of length {} (k = {})",
                node_index,
                sequence,
                sequence.len(),
                k
            );

            let pre_plus: AsciiVectorGenome = sequence.prefix(k - 1).iter().copied().collect();
            let pre_minus: AsciiVectorGenome =
                sequence.suffix(k - 1).reverse_complement_iter().collect();
            let succ_plus: AsciiVectorGenome = sequence.suffix(k - 1).iter().copied().collect();
            let succ_minus: AsciiVectorGenome =
                sequence.prefix(k - 1).reverse_complement_iter().collect();

            let pre_plus = get_or_create_node(&mut bigraph, &mut id_map, pre_plus);
            let pre_minus = get_or_create_node(&mut bigraph, &mut id_map, pre_minus);
            let succ_plus = get_or_create_node(&mut bigraph, &mut id_map, succ_plus);
            let succ_minus = get_or_create_node(&mut bigraph, &mut id_map, succ_minus);

            //println!("Adding edge ({}, {}) and reverse ({}, {})", pre_plus.as_usize(), succ_plus.as_usize(), pre_minus.as_usize(), succ_minus.as_usize());
            bigraph.add_edge(pre_plus, succ_plus, edge_data);
            bigraph.add_edge(pre_minus, succ_minus, reverse_edge_data);
        } else if line.starts_with('L') {
            assert_ne!(k, usize::max_value());

            // Since we are using a hashtable to find the nodes, we can ignore the edges.
        }
    }

    //println!("{:?}", bigraph);
    assert!(header.is_some(), "GFA file has no header");
    assert!(bigraph.verify_node_pairing());
    assert!(bigraph.verify_edge_mirror_property());
    Ok((bigraph, k, header.unwrap()))
}

#[cfg(test)]
mod tests {
    use crate::io::gfa::{read_gfa_as_edge_centric_bigraph, PetGfaEdgeGraph};
    use std::io::BufReader;

    #[test]
    fn test_read_gfa_as_edge_centric_bigraph_simple() {
        let gfa = "H\tKL:Z:3\nS\t1\tACGA\nS\t2\tTCGT";
        let (_bigraph, k, _gfa_header): (PetGfaEdgeGraph<(), ()>, _, _) =
            read_gfa_as_edge_centric_bigraph(BufReader::new(gfa.as_bytes()), false).unwrap();
        assert_eq!(k, 3);
    }
}
