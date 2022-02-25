use crate::io::gfa::BidirectedGfaEdgeData;
use crate::io::SequenceData;
use bigraph::interface::dynamic_bigraph::{DynamicBigraph, DynamicEdgeCentricBigraph};
use bigraph::interface::BidirectedData;
use bigraph::traitgraph::index::GraphIndex;
use bigraph::traitgraph::interface::{GraphBase, ImmutableGraphContainer, StaticGraph};
use bigraph::traitgraph::traitsequence::interface::Sequence;
use bigraph::traitgraph::walks::{EdgeWalk, NodeWalk};
use bio::io::fasta::Record;
use compact_genome::implementation::bit_vec_sequence::BitVectorGenome;
use compact_genome::implementation::DefaultGenome;
use compact_genome::interface::alphabet::Alphabet;
use compact_genome::interface::sequence::{
    EditableGenomeSequence, GenomeSequence, OwnedGenomeSequence,
};
use compact_genome::interface::sequence_store::SequenceStore;
use std::collections::HashMap;
use std::fmt::{Debug, Write};
use std::hash::Hash;
use std::path::Path;

error_chain! {
    foreign_links {
        // For some weird reasons I don't understand, the doc comments have to be put after the item in this macro...
        Io(std::io::Error)
        /// An IO error.
        ;
        Fmt(std::fmt::Error)
        /// An error encountered while trying to format a structure as string.
        ;
        Anyhow(anyhow::Error)
        /// An error passed through anyhow.
        ;
    }

    errors {
        /// A walk is empty.
        EmptyWalkError {
            description("walk is empty")
            display("walk is empty")
        }

        /// An edge has no mirror.
        EdgeWithoutMirror {
            description("an edge has no mirror")
            display("an edge has no mirror")
        }
    }
}

/// Data that can be output as fasta record.
pub trait FastaData<AlphabetType: Alphabet, SourceSequenceStore: SequenceStore<AlphabetType>> {
    /// The type storing the genome sequence of this fasta record.
    type Genome: for<'a> OwnedGenomeSequence<'a, AlphabetType, Self::GenomeSubsequence>
        + for<'a> EditableGenomeSequence<'a, AlphabetType, Self::GenomeSubsequence>;
    /// The subsequence type of `Genome`.
    type GenomeSubsequence: for<'a> GenomeSequence<'a, AlphabetType, Self::GenomeSubsequence>
        + ?Sized;

    /// Returns the sequence of this fasta record.
    fn sequence<'a>(
        &self,
        source_sequence_store: &'a SourceSequenceStore,
    ) -> &'a Self::GenomeSubsequence;
}

/// Write a sequence of walks in a graph as fasta records.
pub fn write_walks_as_fasta<
    'ws,
    AlphabetType: Alphabet + 'static,
    SourceSequenceStore: SequenceStore<AlphabetType>,
    EdgeData: SequenceData<AlphabetType, SourceSequenceStore>,
    Graph: ImmutableGraphContainer<EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    Writer: std::io::Write,
>(
    graph: &Graph,
    source_sequence_store: &SourceSequenceStore,
    kmer_size: usize,
    walks: WalkSource,
    writer: &mut bio::io::fasta::Writer<Writer>,
) -> crate::error::Result<()> {
    for (i, walk) in walks.into_iter().enumerate() {
        if walk.is_empty() {
            return Err(Error::from_kind(ErrorKind::EmptyWalkError).into());
        }

        let mut sequence: DefaultGenome<AlphabetType> = graph
            .edge_data(walk[0])
            .sequence_owned(source_sequence_store);
        for edge in walk.iter().skip(1) {
            let edge_data = graph.edge_data(*edge);
            if let Some(sequence_ref) = edge_data.sequence_ref(source_sequence_store) {
                let sequence_ref = sequence_ref.iter().skip(kmer_size - 1);
                sequence.extend(sequence_ref.cloned());
            } else {
                let sequence_owned: DefaultGenome<AlphabetType> =
                    edge_data.sequence_owned(source_sequence_store);
                let sequence_owned = sequence_owned.iter().skip(kmer_size - 1);
                sequence.extend(sequence_owned.cloned());
            }
        }

        let record =
            bio::io::fasta::Record::with_attrs(&format!("{}", i), None, &sequence.clone_as_vec());
        writer.write_record(&record).map_err(Error::from)?;
    }

    Ok(())
}

/// Write a sequence of walks in a graph as fasta records to a file.
/// The given file is created if it does not exist or truncated if it does exist.
pub fn write_walks_as_fasta_file<
    'ws,
    AlphabetType: Alphabet + 'static,
    SourceSequenceStore: SequenceStore<AlphabetType>,
    EdgeData: SequenceData<AlphabetType, SourceSequenceStore>,
    Graph: ImmutableGraphContainer<EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    P: AsRef<Path>,
>(
    graph: &Graph,
    source_sequence_store: &SourceSequenceStore,
    kmer_size: usize,
    walks: WalkSource,
    path: P,
) -> crate::error::Result<()> {
    write_walks_as_fasta(
        graph,
        source_sequence_store,
        kmer_size,
        walks,
        &mut bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

/// Write a sequence of node-centric walks in a graph as fasta records.
pub fn write_node_centric_walks_as_fasta<
    'ws,
    AlphabetType: Alphabet + 'static,
    SourceSequenceStore: SequenceStore<AlphabetType>,
    NodeData: SequenceData<AlphabetType, SourceSequenceStore>,
    Graph: ImmutableGraphContainer<NodeData = NodeData>,
    Walk: 'ws + for<'w> NodeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> NodeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    Writer: std::io::Write,
>(
    graph: &Graph,
    source_sequence_store: &SourceSequenceStore,
    kmer_size: usize,
    walks: WalkSource,
    writer: &mut bio::io::fasta::Writer<Writer>,
) -> crate::error::Result<()> {
    for (i, walk) in walks.into_iter().enumerate() {
        if walk.is_empty() {
            return Err(Error::from_kind(ErrorKind::EmptyWalkError).into());
        }

        let mut sequence: DefaultGenome<AlphabetType> = graph
            .node_data(walk[0])
            .sequence_owned(source_sequence_store);
        for node in walk.iter().skip(1) {
            let node_data = graph.node_data(*node);
            if let Some(sequence_ref) = node_data.sequence_ref(source_sequence_store) {
                let sequence_ref = sequence_ref.iter().skip(kmer_size - 1);
                sequence.extend(sequence_ref.cloned());
            } else {
                let sequence_owned: DefaultGenome<AlphabetType> =
                    node_data.sequence_owned(source_sequence_store);
                let sequence_owned = sequence_owned.iter().skip(kmer_size - 1);
                sequence.extend(sequence_owned.cloned());
            }
        }

        let record =
            bio::io::fasta::Record::with_attrs(&format!("{}", i), None, &sequence.clone_as_vec());
        writer.write_record(&record).map_err(Error::from)?;
    }

    Ok(())
}

/// Write a sequence of node-centric walks in a graph as fasta records to a file.
/// The given file is created if it does not exist or truncated if it does exist.
pub fn write_node_centric_walks_as_fasta_file<
    'ws,
    AlphabetType: Alphabet + 'static,
    SourceSequenceStore: SequenceStore<AlphabetType>,
    NodeData: SequenceData<AlphabetType, SourceSequenceStore>,
    Graph: ImmutableGraphContainer<NodeData = NodeData>,
    Walk: 'ws + for<'w> NodeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> NodeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    P: AsRef<Path>,
>(
    graph: &Graph,
    source_sequence_store: &SourceSequenceStore,
    kmer_size: usize,
    walks: WalkSource,
    path: P,
) -> crate::error::Result<()> {
    write_node_centric_walks_as_fasta(
        graph,
        source_sequence_store,
        kmer_size,
        walks,
        &mut bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

/// Write a sequence of node-centric walks in a graph as fasta records.
/// The overlaps between the nodes are given by the edges.
pub fn write_node_centric_walks_with_variable_overlaps_as_fasta<
    'ws,
    AlphabetType: Alphabet + 'static,
    SourceSequenceStore: SequenceStore<AlphabetType>,
    NodeData: SequenceData<AlphabetType, SourceSequenceStore>,
    Graph: StaticGraph<NodeData = NodeData, EdgeData = BidirectedGfaEdgeData<()>>,
    Walk: 'ws + for<'w> NodeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> NodeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    Writer: std::io::Write,
>(
    graph: &Graph,
    source_sequence_store: &SourceSequenceStore,
    walks: WalkSource,
    writer: &mut bio::io::fasta::Writer<Writer>,
) -> crate::error::Result<()> {
    for (i, walk) in walks.into_iter().enumerate() {
        if walk.is_empty() {
            return Err(Error::from_kind(ErrorKind::EmptyWalkError).into());
        }

        let mut sequence: DefaultGenome<AlphabetType> = graph
            .node_data(walk[0])
            .sequence_owned(source_sequence_store);
        for (previous_node, node) in walk.iter().take(walk.len() - 1).zip(walk.iter().skip(1)) {
            let node_data = graph.node_data(*node);
            let edge = graph.edges_between(*previous_node, *node).next().unwrap();
            let edge_data = graph.edge_data(edge);
            if let Some(sequence_ref) = node_data.sequence_ref(source_sequence_store) {
                let sequence_ref = sequence_ref.iter().skip(edge_data.overlap);
                sequence.extend(sequence_ref.cloned());
            } else {
                let sequence_owned: DefaultGenome<AlphabetType> =
                    node_data.sequence_owned(source_sequence_store);
                let sequence_owned = sequence_owned.iter().skip(edge_data.overlap);
                sequence.extend(sequence_owned.cloned());
            }
        }

        let record =
            bio::io::fasta::Record::with_attrs(&format!("{}", i), None, &sequence.clone_as_vec());
        writer.write_record(&record).map_err(Error::from)?;
    }

    Ok(())
}

/// Write a sequence of node-centric walks in a graph as fasta records to a file.
/// The overlaps between the nodes are given by the edges.
/// The given file is created if it does not exist or truncated if it does exist.
pub fn write_node_centric_walks_with_variable_overlaps_as_fasta_file<
    'ws,
    AlphabetType: Alphabet + 'static,
    SourceSequenceStore: SequenceStore<AlphabetType>,
    NodeData: SequenceData<AlphabetType, SourceSequenceStore>,
    Graph: StaticGraph<NodeData = NodeData, EdgeData = BidirectedGfaEdgeData<()>>,
    Walk: 'ws + for<'w> NodeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> NodeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    P: AsRef<Path>,
>(
    graph: &Graph,
    source_sequence_store: &SourceSequenceStore,
    walks: WalkSource,
    path: P,
) -> crate::error::Result<()> {
    write_node_centric_walks_with_variable_overlaps_as_fasta(
        graph,
        source_sequence_store,
        walks,
        &mut bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

/////////////////////////////
////// NODE CENTRIC IO //////
/////////////////////////////

/// Raw data of a fasta record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastaNodeData<SequenceHandle> {
    /// The id of the fasta record.
    pub id: String,
    /// The optional description of the fasta record.
    pub description: Option<String>,
    /// A handle of the sequence of the fasta record.
    pub sequence_handle: SequenceHandle,
    /// True if this is the forwards variant of the fasta record, false if it is the backwards variant.
    /// This is always true for records read from a fasta file, but the mirror nodes of the fasta records have this set to false.
    pub forwards: bool,
}

impl<SequenceHandle: Clone> BidirectedData for FastaNodeData<SequenceHandle> {
    fn mirror(&self) -> Self {
        let mut result = self.clone();
        result.forwards = !result.forwards;
        result
    }
}

fn parse_fasta_record<AlphabetType: Alphabet, GenomeSequenceStore: SequenceStore<AlphabetType>>(
    record: Record,
    target_sequence_store: &mut GenomeSequenceStore,
) -> Result<FastaNodeData<GenomeSequenceStore::Handle>> {
    let id = record.id().to_owned();
    let description = record.desc().map(ToOwned::to_owned);
    let sequence_handle = target_sequence_store
        .add_from_slice_u8(record.seq())
        .unwrap_or_else(|error| panic!("Genome sequence with id {id} is invalid: {error:?}"));
    Ok(FastaNodeData {
        id,
        description,
        sequence_handle,
        forwards: true,
    })
}

/*
/// Read a genome graph in fasta format into a node-centric representation from a file.
pub fn read_bigraph_from_fasta_as_node_centric_from_file<
    P: AsRef<Path>,
    GenomeSequenceStore: SequenceStore,
    NodeData: From<FastaNodeData<GenomeSequenceStore::Handle>> + BidirectedData,
    EdgeData: Default + Clone,
    Graph: DynamicNodeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    path: P,
    target_sequence_store: &mut GenomeSequenceStore,
) -> crate::error::Result<Graph> {
    read_bigraph_from_fasta_as_node_centric(
        bio::io::fasta::Reader::from_file(path).map_err(Error::from)?,
        target_sequence_store,
    )
}

/// Read a genome graph in fasta format into a node-centric representation.
pub fn read_bigraph_from_fasta_as_node_centric<
    R: std::io::Read,
    GenomeSequenceStore: SequenceStore,
    NodeData: From<FastaNodeData<GenomeSequenceStore::Handle>> + BidirectedData,
    EdgeData: Default + Clone,
    Graph: DynamicNodeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    reader: bio::io::fasta::Reader<R>,
    target_sequence_store: &mut GenomeSequenceStore,
) -> crate::error::Result<Graph> {
    let mut bigraph = Graph::default();
    let mut edges = Vec::new();

    for record in reader.records() {
        let record =
            parse_fasta_record(record.map_err(Error::from)?, target_sequence_store)?;
        edges.extend(record.edges.iter().map(|e| BiEdge {
            from_node: record.id,
            plain_edge: e.clone(),
        }));
        let record_id = record.id;
        let id = bigraph.add_node(record.into());
        debug_assert_eq!(id, record_id.into());
    }

    bigraph.add_mirror_nodes();
    debug_assert!(bigraph.verify_node_pairing());

    for edge in edges {
        let from_node = if edge.plain_edge.from_side {
            edge.from_node.into()
        } else {
            bigraph.mirror_node(edge.from_node.into()).unwrap()
        };
        let to_node = if edge.plain_edge.to_side {
            edge.plain_edge.to_node.into()
        } else {
            bigraph.mirror_node(edge.plain_edge.to_node.into()).unwrap()
        };
        bigraph.add_edge(from_node, to_node, EdgeData::default());
    }

    bigraph.add_node_centric_mirror_edges();
    debug_assert!(bigraph.verify_node_mirror_property());
    Ok(bigraph)
}

/// Write a genome graph in bcalm2 fasta format from a node-centric representation to a file.
pub fn write_node_centric_bigraph_to_bcalm2_to_file<
    P: AsRef<Path>,
    GenomeSequenceStore: SequenceStore,
    NodeData, //: Into<PlainBCalm2NodeData<IndexType>>,
    EdgeData: Default + Clone,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    graph: &Graph,
    source_sequence_store: &GenomeSequenceStore,
    path: P,
) -> crate::error::Result<()>
    where
            for<'a> PlainBCalm2NodeData<GenomeSequenceStore::Handle>: From<&'a NodeData>,
{
    write_node_centric_bigraph_to_bcalm2(
        graph,
        source_sequence_store,
        bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

/// Write a genome graph in bcalm2 fasta format from a node-centric representation.
pub fn write_node_centric_bigraph_to_bcalm2<
    W: std::io::Write,
    GenomeSequenceStore: SequenceStore,
    NodeData,
    EdgeData: Default + Clone,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    graph: &Graph,
    source_sequence_store: &GenomeSequenceStore,
    mut writer: bio::io::fasta::Writer<W>,
) -> crate::error::Result<()>
    where
            for<'a> PlainBCalm2NodeData<GenomeSequenceStore::Handle>: From<&'a NodeData>,
{
    let mut output_nodes = vec![false; graph.node_count()];

    for node_id in graph.node_indices() {
        if !output_nodes[graph
            .mirror_node(node_id)
            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutMirror))?
            .as_usize()]
        {
            output_nodes[node_id.as_usize()] = true;
        }
    }

    for node_id in graph.node_indices() {
        if output_nodes[node_id.as_usize()] {
            let node_data = PlainBCalm2NodeData::from(graph.node_data(node_id));
            let mirror_node_id = graph
                .mirror_node(node_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutMirror))?;
            /*let mirror_node_data = PlainBCalm2NodeData::<IndexType>::from(
                graph
                    .node_data(mirror_node_id)
                    .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutMirror))?,
            );*/
            let mut out_neighbors_plus = Vec::new();
            let mut out_neighbors_minus = Vec::new();

            for neighbor in graph.out_neighbors(node_id) {
                let neighbor_node_id = neighbor.node_id.as_usize();

                out_neighbors_plus.push((
                    true,
                    if output_nodes[neighbor_node_id] {
                        neighbor.node_id.as_usize()
                    } else {
                        graph
                            .mirror_node(neighbor.node_id)
                            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutMirror))?
                            .as_usize()
                    },
                    output_nodes[neighbor_node_id],
                ));
            }
            for neighbor in graph.out_neighbors(mirror_node_id) {
                let neighbor_node_id = neighbor.node_id.as_usize();

                out_neighbors_minus.push((
                    false,
                    if output_nodes[neighbor_node_id] {
                        neighbor.node_id.as_usize()
                    } else {
                        graph
                            .mirror_node(neighbor.node_id)
                            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutMirror))?
                            .as_usize()
                    },
                    output_nodes[neighbor_node_id],
                ));
            }

            out_neighbors_plus.sort_unstable();
            out_neighbors_minus.sort_unstable();
            out_neighbors_plus.append(&mut out_neighbors_minus);
            let out_neighbors = out_neighbors_plus;

            let mut printed_node_id = String::new();
            write!(printed_node_id, "{}", node_data.id).map_err(Error::from)?;
            let node_description =
                write_plain_bcalm2_node_data_to_bcalm2(&node_data, out_neighbors)?;
            let node_sequence = source_sequence_store
                .get(&node_data.sequence_handle)
                .clone_as_vec();

            writer
                .write(&printed_node_id, Some(&node_description), &node_sequence)
                .map_err(Error::from)?;
        }
    }

    Ok(())
}
*/

/////////////////////////////
////// EDGE CENTRIC IO //////
/////////////////////////////

/// Read a genome graph in fasta format into an edge-centric representation from a file.
pub fn read_bigraph_from_fasta_as_edge_centric_from_file<
    P: AsRef<Path> + Debug,
    AlphabetType: Alphabet + Hash + Eq + Clone + 'static,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
    NodeData: Default + Clone,
    EdgeData: From<FastaNodeData<GenomeSequenceStore::Handle>> + Clone + Eq + BidirectedData,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    path: P,
    target_sequence_store: &mut GenomeSequenceStore,
    kmer_size: usize,
) -> crate::error::Result<Graph>
where
    <GenomeSequenceStore as SequenceStore<AlphabetType>>::Handle: Clone,
{
    read_bigraph_from_fasta_as_edge_centric(
        bio::io::fasta::Reader::from_file(path).map_err(Error::from)?,
        target_sequence_store,
        kmer_size,
    )
}

fn get_or_create_node<
    Graph: DynamicBigraph,
    AlphabetType: Alphabet,
    Genome: for<'a> OwnedGenomeSequence<'a, AlphabetType, GenomeSubsequence> + Hash + Eq + Clone,
    GenomeSubsequence: for<'a> GenomeSequence<'a, AlphabetType, GenomeSubsequence> + ?Sized,
>(
    bigraph: &mut Graph,
    id_map: &mut HashMap<Genome, <Graph as GraphBase>::NodeIndex>,
    genome: Genome,
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

/// Read a genome graph in fasta format into an edge-centric representation.
pub fn read_bigraph_from_fasta_as_edge_centric<
    R: std::io::BufRead,
    AlphabetType: Alphabet + Hash + Eq + Clone + 'static,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
    NodeData: Default + Clone,
    EdgeData: From<FastaNodeData<GenomeSequenceStore::Handle>> + Clone + Eq + BidirectedData,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    reader: bio::io::fasta::Reader<R>,
    target_sequence_store: &mut GenomeSequenceStore,
    kmer_size: usize,
) -> crate::error::Result<Graph>
where
    <Graph as GraphBase>::NodeIndex: Clone,
    <GenomeSequenceStore as SequenceStore<AlphabetType>>::Handle: Clone,
{
    let mut bigraph = Graph::default();
    let mut id_map = HashMap::new();
    let node_kmer_size = kmer_size - 1;

    for record in reader.records() {
        let record: FastaNodeData<GenomeSequenceStore::Handle> =
            parse_fasta_record(record.map_err(Error::from)?, target_sequence_store)?;
        let sequence = target_sequence_store.get(&record.sequence_handle);
        let prefix = sequence.prefix(node_kmer_size);
        let suffix = sequence.suffix(node_kmer_size);

        let pre_plus: BitVectorGenome<AlphabetType> = prefix.convert();
        let pre_minus: BitVectorGenome<AlphabetType> = suffix.convert_with_reverse_complement();
        let succ_plus: BitVectorGenome<AlphabetType> = suffix.convert();
        let succ_minus: BitVectorGenome<AlphabetType> = prefix.convert_with_reverse_complement();

        let pre_plus = get_or_create_node(&mut bigraph, &mut id_map, pre_plus);
        let pre_minus = get_or_create_node(&mut bigraph, &mut id_map, pre_minus);
        let succ_plus = get_or_create_node(&mut bigraph, &mut id_map, succ_plus);
        let succ_minus = get_or_create_node(&mut bigraph, &mut id_map, succ_minus);

        bigraph.add_edge(pre_plus, succ_plus, record.clone().into());
        bigraph.add_edge(pre_minus, succ_minus, record.mirror().into());
    }

    debug_assert!(bigraph.verify_node_pairing());
    debug_assert!(bigraph.verify_edge_mirror_property());
    Ok(bigraph)
}

/// Write a genome graph in fasta format from an edge-centric representation to a file.
pub fn write_edge_centric_bigraph_to_fasta_to_file<
    P: AsRef<Path>,
    AlphabetType: Alphabet,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
    NodeData, //: Into<PlainBCalm2NodeData<IndexType>>,
    EdgeData: BidirectedData + Clone + Eq,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    graph: &Graph,
    source_sequence_store: &GenomeSequenceStore,
    path: P,
) -> crate::error::Result<()>
where
    for<'a> FastaNodeData<GenomeSequenceStore::Handle>: From<&'a EdgeData>,
{
    write_edge_centric_bigraph_to_fasta(
        graph,
        source_sequence_store,
        bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

/// Write a genome graph in fasta format from an edge-centric representation.
pub fn write_edge_centric_bigraph_to_fasta<
    W: std::io::Write,
    AlphabetType: Alphabet,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
    NodeData,
    EdgeData: BidirectedData + Clone + Eq,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    graph: &Graph,
    source_sequence_store: &GenomeSequenceStore,
    mut writer: bio::io::fasta::Writer<W>,
) -> crate::error::Result<()>
where
    for<'a> FastaNodeData<GenomeSequenceStore::Handle>: From<&'a EdgeData>,
{
    let mut output_edges = vec![false; graph.edge_count()];

    for edge_id in graph.edge_indices() {
        if !output_edges[graph
            .mirror_edge_edge_centric(edge_id)
            .ok_or_else(|| Error::from(ErrorKind::EdgeWithoutMirror))?
            .as_usize()]
        {
            output_edges[edge_id.as_usize()] = true;
        }
    }

    for edge_id in graph.edge_indices() {
        if output_edges[edge_id.as_usize()] {
            let node_data = FastaNodeData::from(graph.edge_data(edge_id));

            let mut printed_node_id = String::new();
            write!(printed_node_id, "{}", node_data.id).map_err(Error::from)?;
            let node_sequence = source_sequence_store
                .get(&node_data.sequence_handle)
                .clone_as_vec();

            writer
                .write(
                    &printed_node_id,
                    node_data.description.as_deref(),
                    &node_sequence,
                )
                .map_err(Error::from)?;
        }
    }

    Ok(())
}

//////////////////////////////////////
////// PARALLEL EDGE CENTRIC IO //////
//////////////////////////////////////

/*
/// Read a genome graph in fasta format into an edge-centric representation from a file.
pub fn read_bigraph_from_fasta_as_edge_centric_from_file_in_parallel<
    P: AsRef<Path> + Debug,
    GenomeSequenceStore: SequenceStore + Send,
    NodeData: Default + Clone,
    EdgeData: From<FastaNodeData<GenomeSequenceStore::Handle>> + Clone + Eq + BidirectedData,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default + Send,
>(
    path: P,
    target_sequence_store: &mut GenomeSequenceStore,
    kmer_size: usize,
) -> crate::error::Result<Graph>
where
    <Graph as GraphBase>::NodeIndex: Clone + Send + Sync,
    <GenomeSequenceStore as SequenceStore>::Handle: Clone,
{
    read_bigraph_from_fasta_as_edge_centric_in_parallel(
        bio::io::fasta::Reader::from_file(path).map_err(Error::from)?,
        target_sequence_store,
        kmer_size,
    )
}

fn get_or_create_node_in_parallel<
    Graph: DynamicBigraph,
    Genome: for<'a> OwnedGenomeSequence<'a, GenomeSubsequence> + Hash + Eq + Clone,
    GenomeSubsequence: for<'a> GenomeSequence<'a, GenomeSubsequence> + ?Sized,
>(
    bigraph: &Mutex<Graph>,
    id_map: &DashMap<Genome, <Graph as GraphBase>::NodeIndex>,
    genome: Genome,
) -> <Graph as GraphBase>::NodeIndex
where
    <Graph as GraphBase>::NodeData: Default,
    <Graph as GraphBase>::EdgeData: Clone,
{
    if let Some(node) = id_map.get(&genome) {
        *node
    } else {
        let bigraph = &mut *bigraph.lock().unwrap();
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

/// Read a genome graph in fasta format into an edge-centric representation.
pub fn read_bigraph_from_fasta_as_edge_centric_in_parallel<
    R: std::io::BufRead + Send,
    GenomeSequenceStore: SequenceStore + Send,
    NodeData: Default + Clone,
    EdgeData: From<FastaNodeData<GenomeSequenceStore::Handle>> + Clone + Eq + BidirectedData,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default + Send,
>(
    reader: bio::io::fasta::Reader<R>,
    target_sequence_store: &mut GenomeSequenceStore,
    kmer_size: usize,
) -> crate::error::Result<Graph>
where
    <Graph as GraphBase>::NodeIndex: Clone + Send + Sync,
    <GenomeSequenceStore as SequenceStore>::Handle: Clone,
{
    todo!("This results in broken graphs");

    let bigraph = Mutex::new(Graph::default());
    let target_sequence_store = Mutex::new(target_sequence_store);
    let id_map = DashMap::new();
    let node_kmer_size = kmer_size - 1;

    let result: Result<()> = reader
        .records()
        .par_bridge()
        .map(|record| {
            let mut locked_target_sequence_store = target_sequence_store.lock().unwrap();
            let record: FastaNodeData<GenomeSequenceStore::Handle> =
                parse_fasta_record(record?, *locked_target_sequence_store)?;
            let sequence = locked_target_sequence_store.get(&record.sequence_handle);
            let prefix = sequence.prefix(node_kmer_size);
            let suffix = sequence.suffix(node_kmer_size);

            let pre_plus: BitVectorGenome = prefix.convert();
            let pre_minus: BitVectorGenome = suffix.convert_with_reverse_complement();
            let succ_plus: BitVectorGenome = suffix.convert();
            let succ_minus: BitVectorGenome = prefix.convert_with_reverse_complement();
            drop(locked_target_sequence_store);

            let pre_plus = get_or_create_node_in_parallel(&bigraph, &id_map, pre_plus);
            let pre_minus = get_or_create_node_in_parallel(&bigraph, &id_map, pre_minus);
            let succ_plus = get_or_create_node_in_parallel(&bigraph, &id_map, succ_plus);
            let succ_minus = get_or_create_node_in_parallel(&bigraph, &id_map, succ_minus);

            let bigraph = &mut *bigraph.lock().unwrap();
            bigraph.add_edge(pre_plus, succ_plus, record.clone().into());
            bigraph.add_edge(pre_minus, succ_minus, record.mirror().into());
            Result::<()>::Ok(())
        })
        .collect();
    result?;

    let bigraph = bigraph.into_inner().unwrap();
    debug_assert!(bigraph.verify_node_pairing());
    debug_assert!(bigraph.verify_edge_mirror_property());
    Ok(bigraph)
}
*/
