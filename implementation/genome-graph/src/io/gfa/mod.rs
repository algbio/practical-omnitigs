use crate::error::{Error, ErrorKind, Result};
use crate::io::SequenceData;
use bigraph::interface::dynamic_bigraph::{DynamicBigraph, DynamicEdgeCentricBigraph};
use bigraph::interface::BidirectedData;
use bigraph::traitgraph::index::GraphIndex;
use bigraph::traitgraph::interface::GraphBase;
use compact_genome::implementation::DefaultGenome;
use compact_genome::interface::alphabet::Alphabet;
use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
use compact_genome::interface::sequence_store::SequenceStore;
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader};
use std::path::Path;
use traitgraph_algo::dijkstra::DijkstraWeightedEdgeData;

/// Type of graphs read from gfa files.
pub type PetGfaGraph<NodeData, EdgeData, SequenceHandle> =
    crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            BidirectedGfaNodeData<SequenceHandle, NodeData>,
            BidirectedGfaEdgeData<EdgeData>,
            usize,
        >,
    >;

/// The edge-centric variant of the type of graphs read from gfa files.
pub type PetGfaEdgeGraph<NodeData, EdgeData, SequenceHandle> =
    crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            NodeData,
            BidirectedGfaNodeData<SequenceHandle, EdgeData>,
            usize,
        >,
    >;

/// Node data of a bidirected graph read from GFA.
#[derive(Eq, PartialEq, Debug, Clone, Default)]
pub struct BidirectedGfaNodeData<SequenceHandle, Data> {
    /// The sequence of this node. If forward is false, then this must be reverse complemented.
    pub sequence_handle: SequenceHandle,
    /// True if this node is the forward node of sequence, false if it is the reverse complement node.
    pub forward: bool,
    /// Further data.
    pub data: Data,
}

impl<SequenceHandle: Clone, Data: BidirectedData> BidirectedData
    for BidirectedGfaNodeData<SequenceHandle, Data>
{
    fn mirror(&self) -> Self {
        Self {
            sequence_handle: self.sequence_handle.clone(),
            forward: !self.forward,
            data: self.data.mirror(),
        }
    }
}

impl<SequenceHandle, Data: DijkstraWeightedEdgeData<usize>> DijkstraWeightedEdgeData<usize>
    for BidirectedGfaNodeData<SequenceHandle, Data>
{
    fn weight(&self) -> usize {
        self.data.weight()
    }
}

impl<AlphabetType: Alphabet, GenomeSequenceStore: SequenceStore<AlphabetType>, Data>
    SequenceData<AlphabetType, GenomeSequenceStore>
    for BidirectedGfaNodeData<GenomeSequenceStore::Handle, Data>
{
    fn sequence_handle(&self) -> &GenomeSequenceStore::Handle {
        &self.sequence_handle
    }

    fn sequence_ref<'a>(
        &self,
        source_sequence_store: &'a GenomeSequenceStore,
    ) -> Option<&'a <GenomeSequenceStore as SequenceStore<AlphabetType>>::SequenceRef> {
        if self.forward {
            let handle =
                <BidirectedGfaNodeData<GenomeSequenceStore::Handle, Data> as SequenceData<
                    AlphabetType,
                    GenomeSequenceStore,
                >>::sequence_handle(self);
            Some(source_sequence_store.get(handle))
        } else {
            None
        }
    }

    fn sequence_owned<
        ResultSequence: for<'a> OwnedGenomeSequence<'a, AlphabetType, ResultSubsequence>,
        ResultSubsequence: for<'a> GenomeSequence<'a, AlphabetType, ResultSubsequence> + ?Sized,
    >(
        &self,
        source_sequence_store: &GenomeSequenceStore,
    ) -> ResultSequence {
        let handle = <BidirectedGfaNodeData<GenomeSequenceStore::Handle, Data> as SequenceData<
            AlphabetType,
            GenomeSequenceStore,
        >>::sequence_handle(self);
        if self.forward {
            source_sequence_store.get(handle).convert()
        } else {
            source_sequence_store
                .get(handle)
                .convert_with_reverse_complement()
        }
    }
}

/// Edge data of a bidirected graph read from GFA.
#[derive(Eq, PartialEq, Debug, Clone, Default)]
pub struct BidirectedGfaEdgeData<Data> {
    /// Size of the overlap between the tail and head nodes.
    pub overlap: usize,
    /// Further data.
    pub data: Data,
}

/// Properties of a GFA file that was read.
pub struct GfaReadFileProperties {
    /// The order of the node-centric de Bruijn graph stored in the GFA file. If the GFA file does not contain the respective header field, then this field is usize::max_value().
    pub k: usize,

    /// The header of the GFA file. Should the GFA file have multiple header lines, it is undefined which line is reported. If the GFA file has no header lines, then this field is None.
    pub header: Option<String>,
}

/// Read a bigraph in gfa format from a file.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_bigraph_from_file<
    P: AsRef<Path>,
    AlphabetType: Alphabet,
    GenomeSequenceStoreHandle: Clone,
    GenomeSequenceStoreRef: for<'a> GenomeSequence<'a, AlphabetType, GenomeSequenceStoreRef> + Debug + ?Sized,
    GenomeSequenceStore: SequenceStore<
        AlphabetType,
        Handle = GenomeSequenceStoreHandle,
        SequenceRef = GenomeSequenceStoreRef,
    >,
    NodeData: From<BidirectedGfaNodeData<GenomeSequenceStore::Handle, ()>>,
    EdgeData: From<BidirectedGfaEdgeData<()>>,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    gfa_file: P,
    target_sequence_store: &mut GenomeSequenceStore,
    ignore_k: bool,
    allow_messy_edges: bool,
) -> Result<(Graph, GfaReadFileProperties)> {
    read_gfa_as_bigraph(
        BufReader::new(File::open(gfa_file)?),
        target_sequence_store,
        ignore_k,
        allow_messy_edges,
    )
}

/// Read a bigraph in gfa format from a `BufRead`.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_bigraph<
    R: BufRead,
    AlphabetType: Alphabet,
    GenomeSequenceStoreHandle: Clone,
    GenomeSequenceStoreRef: for<'a> GenomeSequence<'a, AlphabetType, GenomeSequenceStoreRef> + Debug + ?Sized,
    GenomeSequenceStore: SequenceStore<
        AlphabetType,
        Handle = GenomeSequenceStoreHandle,
        SequenceRef = GenomeSequenceStoreRef,
    >,
    NodeData: From<BidirectedGfaNodeData<GenomeSequenceStore::Handle, ()>>,
    EdgeData: From<BidirectedGfaEdgeData<()>>,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    gfa: R,
    target_sequence_store: &mut GenomeSequenceStore,
    ignore_k: bool,
    allow_messy_edges: bool,
) -> Result<(Graph, GfaReadFileProperties)> {
    let mut graph = Graph::default();
    let mut k = usize::max_value();
    let mut header = None;
    let mut node_name_map = HashMap::new();

    for line in gfa.lines() {
        let line = line?;

        if line.starts_with('H') {
            debug_assert!(graph.is_empty());
            header = Some(line.to_owned());
            for column in line.split('\t') {
                if let Some(stripped) = column.strip_prefix("KL:Z:") {
                    debug_assert_eq!(k, usize::max_value());
                    k = stripped.parse().unwrap();
                }
            }
        } else if line.starts_with('S') {
            if !allow_messy_edges {
                debug_assert_eq!(graph.edge_count(), 0);
            }
            if !ignore_k {
                debug_assert_ne!(k, usize::max_value());
            }

            let mut columns = line.split('\t').skip(1);
            let node_name: &str = columns.next().unwrap();

            let sequence = columns.next().unwrap().as_bytes();
            let sequence_handle = target_sequence_store
                .add_from_slice_u8(sequence)
                .unwrap_or_else(|error| {
                    panic!("Genome sequence with node_name {node_name} is invalid: {error:?}")
                });
            let sequence = target_sequence_store.get(&sequence_handle);
            debug_assert!(
                sequence.len() >= k || ignore_k,
                "Node {} has sequence '{:?}' of length {} (k = {})",
                node_name,
                sequence,
                sequence.len(),
                k
            );

            let n1 = graph.add_node(
                BidirectedGfaNodeData {
                    sequence_handle: sequence_handle.clone(),
                    forward: true,
                    data: Default::default(),
                }
                .into(),
            );
            let n2 = graph.add_node(
                BidirectedGfaNodeData {
                    sequence_handle: sequence_handle.clone(),
                    forward: false,
                    data: Default::default(),
                }
                .into(),
            );
            graph.set_mirror_nodes(n1, n2);
            node_name_map.insert(node_name.to_owned(), n1);
        } else if line.starts_with('L') {
            if !ignore_k {
                debug_assert_ne!(k, usize::max_value());
            }

            let mut columns = line.split('\t').skip(1);
            let n1_name = columns.next().unwrap();
            let n1_direction = if columns.next().unwrap() == "+" { 0 } else { 1 };
            let n2_name = columns.next().unwrap();
            let n2_direction = if columns.next().unwrap() == "+" { 0 } else { 1 };
            let overlap = if let Some(overlap) = columns.next() {
                if let Some(overlap) = overlap.strip_suffix('M') {
                    overlap
                        .parse()
                        .map_err(|_| Error::from_kind(ErrorKind::GfaUnknownOverlapPattern))?
                } else {
                    return Err(ErrorKind::GfaUnknownOverlapPattern.into());
                }
            } else {
                return Err(ErrorKind::GfaMissingOverlapPattern.into());
            };

            if let (Some(n1), Some(n2)) = (node_name_map.get(n1_name), node_name_map.get(n2_name)) {
                let n1 = (n1.as_usize() + n1_direction).into();
                let n2 = (n2.as_usize() + n2_direction).into();

                let has_edge = graph.contains_edge_between(n1, n2);
                debug_assert_eq!(
                    has_edge,
                    graph.contains_edge_between(
                        graph.mirror_node(n2).unwrap(),
                        graph.mirror_node(n1).unwrap()
                    )
                );

                if !has_edge {
                    let edge_data = BidirectedGfaEdgeData { data: (), overlap };
                    graph.add_edge(n1, n2, edge_data.clone().into());
                    graph.add_edge(
                        graph.mirror_node(n2).unwrap(),
                        graph.mirror_node(n1).unwrap(),
                        edge_data.into(),
                    );
                }
            } else {
                return Err(ErrorKind::GfaMissingNode.into());
            }
        }
    }

    if ignore_k {
        k = 0;
    }

    Ok((graph, GfaReadFileProperties { k, header }))
}

/// Read an edge-centric bigraph in gfa format from a file.
/// This method also returns the k-mer length given in the gfa file as well as the full gfa header.
pub fn read_gfa_as_edge_centric_bigraph_from_file<
    P: AsRef<Path>,
    AlphabetType: Alphabet + Clone + Eq + Hash + 'static,
    GenomeSequenceStoreHandle: Clone + Eq,
    GenomeSequenceStoreRef: for<'a> GenomeSequence<'a, AlphabetType, GenomeSequenceStoreRef> + Debug + ?Sized,
    GenomeSequenceStore: SequenceStore<
        AlphabetType,
        Handle = GenomeSequenceStoreHandle,
        SequenceRef = GenomeSequenceStoreRef,
    >,
    NodeData: Default,
    EdgeData: Default
        + BidirectedData
        + Eq
        + Clone
        + From<BidirectedGfaNodeData<GenomeSequenceStore::Handle, ()>>,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default + std::fmt::Debug,
>(
    gfa_file: P,
    target_sequence_store: &mut GenomeSequenceStore,
    estimate_k: bool,
) -> Result<(Graph, GfaReadFileProperties)> {
    read_gfa_as_edge_centric_bigraph(
        BufReader::new(File::open(gfa_file)?),
        target_sequence_store,
        estimate_k,
    )
}

fn get_or_create_node<
    Graph: DynamicBigraph,
    AlphabetType: Alphabet,
    G: for<'a> OwnedGenomeSequence<'a, AlphabetType, GenomeSubsequence> + Hash + Eq + Clone,
    GenomeSubsequence: for<'a> GenomeSequence<'a, AlphabetType, GenomeSubsequence> + ?Sized,
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

        let reverse_complement = genome.clone_as_reverse_complement();
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
    AlphabetType: Alphabet + Clone + Eq + Hash + 'static,
    GenomeSequenceStoreHandle: Clone + Eq,
    GenomeSequenceStoreRef: for<'a> GenomeSequence<'a, AlphabetType, GenomeSequenceStoreRef> + Debug + ?Sized,
    GenomeSequenceStore: SequenceStore<
        AlphabetType,
        Handle = GenomeSequenceStoreHandle,
        SequenceRef = GenomeSequenceStoreRef,
    >,
    NodeData: Default,
    EdgeData: Default
        + BidirectedData
        + Eq
        + Clone
        + From<BidirectedGfaNodeData<GenomeSequenceStore::Handle, ()>>,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default + std::fmt::Debug,
>(
    gfa: R,
    target_sequence_store: &mut GenomeSequenceStore,
    estimate_k: bool,
) -> Result<(Graph, GfaReadFileProperties)> {
    debug_assert!(!estimate_k, "Estimating k not supported yet");

    let mut bigraph = Graph::default();
    let mut id_map = HashMap::new();
    let mut k = usize::max_value();
    let mut header = None;

    for line in gfa.lines() {
        let line = line?;

        if line.starts_with('H') {
            debug_assert!(bigraph.is_empty());
            header = Some(line.clone());
            for column in line.split('\t') {
                if let Some(stripped) = column.strip_prefix("KL:Z:") {
                    debug_assert_eq!(k, usize::max_value());
                    k = stripped.parse().unwrap();
                }
            }
        } else if line.starts_with('S') {
            debug_assert_ne!(k, usize::max_value());

            let mut columns = line.split('\t').skip(1);
            let node_index: usize = columns.next().unwrap().parse().unwrap();
            debug_assert_eq!((node_index - 1) * 2, bigraph.edge_count());

            let sequence = columns.next().unwrap().as_bytes();
            //println!("sequence {}", sequence);
            let sequence_handle = target_sequence_store
                .add_from_slice_u8(sequence)
                .unwrap_or_else(|error| {
                    panic!("Genome sequence with node_index {node_index} is invalid: {error:?}")
                });
            let sequence = target_sequence_store.get(&sequence_handle);
            let edge_data = BidirectedGfaNodeData {
                sequence_handle: sequence_handle.clone(),
                forward: true,
                data: Default::default(),
            };
            let edge_data: EdgeData = edge_data.into();
            let reverse_edge_data = edge_data.mirror();

            debug_assert!(columns.next().is_none());
            debug_assert!(
                sequence.len() >= k,
                "Node {} has sequence '{:?}' of length {} (k = {})",
                node_index,
                sequence,
                sequence.len(),
                k
            );

            let pre_plus: DefaultGenome<AlphabetType> = sequence.prefix(k - 1).convert();
            let pre_minus: DefaultGenome<AlphabetType> =
                sequence.suffix(k - 1).reverse_complement_iter().collect();
            let succ_plus: DefaultGenome<AlphabetType> = sequence.suffix(k - 1).convert();
            let succ_minus: DefaultGenome<AlphabetType> =
                sequence.prefix(k - 1).reverse_complement_iter().collect();

            let pre_plus = get_or_create_node(&mut bigraph, &mut id_map, pre_plus);
            let pre_minus = get_or_create_node(&mut bigraph, &mut id_map, pre_minus);
            let succ_plus = get_or_create_node(&mut bigraph, &mut id_map, succ_plus);
            let succ_minus = get_or_create_node(&mut bigraph, &mut id_map, succ_minus);

            //println!("Adding edge ({}, {}) and reverse ({}, {})", pre_plus.as_usize(), succ_plus.as_usize(), pre_minus.as_usize(), succ_minus.as_usize());
            bigraph.add_edge(pre_plus, succ_plus, edge_data);
            bigraph.add_edge(pre_minus, succ_minus, reverse_edge_data);
        } else if line.starts_with('L') {
            debug_assert_ne!(k, usize::max_value());

            // Since we are using a hashtable to find the nodes, we can ignore the edges.
        }
    }

    //println!("{:?}", bigraph);
    debug_assert!(header.is_some(), "GFA file has no header");
    debug_assert!(bigraph.verify_node_pairing());
    debug_assert!(bigraph.verify_edge_mirror_property());
    Ok((bigraph, GfaReadFileProperties { k, header }))
}

#[cfg(test)]
mod tests {
    use crate::io::gfa::{
        read_gfa_as_edge_centric_bigraph, GfaReadFileProperties, PetGfaEdgeGraph,
    };
    use compact_genome::implementation::DefaultSequenceStore;
    use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
    use std::io::BufReader;

    #[test]
    fn test_read_gfa_as_edge_centric_bigraph_simple() {
        let gfa = "H\tKL:Z:3\nS\t1\tACGA\nS\t2\tTCGT";
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();
        let (_bigraph, GfaReadFileProperties { k, .. }): (PetGfaEdgeGraph<(), (), _>, _) =
            read_gfa_as_edge_centric_bigraph(
                BufReader::new(gfa.as_bytes()),
                &mut sequence_store,
                false,
            )
            .unwrap();
        debug_assert_eq!(k, 3);
    }
}
