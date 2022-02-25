use crate::bigraph::interface::dynamic_bigraph::DynamicEdgeCentricBigraph;
use crate::bigraph::interface::dynamic_bigraph::DynamicNodeCentricBigraph;
use crate::io::SequenceData;
use bigraph::interface::{dynamic_bigraph::DynamicBigraph, BidirectedData};
use bigraph::traitgraph::index::GraphIndex;
use bigraph::traitgraph::interface::GraphBase;
use bigraph::traitgraph::traitsequence::interface::Sequence;
use bio::io::fasta::Record;
use compact_genome::implementation::bit_vec_sequence::BitVectorGenome;
use compact_genome::interface::alphabet::Alphabet;
use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
use compact_genome::interface::sequence_store::SequenceStore;
use num_traits::NumCast;
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
        /// Any error passed through anyhow.
        ;
    }

    errors {
        /// A node id that cannot be parsed into `usize`.
        BCalm2IdError(id: String) {
            description("invalid node id")
            display("invalid node id: '{:?}'", id)
        }

        /// The length of a sequence does not match the length given as fasta parameter.
        BCalm2LengthError(length: usize, sequence_length: usize) {
            description("the length in the description of a node does not match the length of its sequence")
            display("the length in the description of a node ({}) does not match the length of its sequence {}", length, sequence_length)
        }

        /// The fasta comment contains an unknown parameter.
        BCalm2UnknownParameterError(parameter: String) {
            description("unknown parameter")
            display("unknown parameter: '{:?}'", parameter)
        }

        /// The fasta comment contains a duplicate parameter.
        BCalm2DuplicateParameterError(parameter: String) {
            description("duplicate parameter")
            display("duplicate parameter: '{:?}'", parameter)
        }

        /// The fasta comment contains a parameter that cannot be parsed.
        BCalm2MalformedParameterError(parameter: String) {
            description("malformed parameter")
            display("malformed parameter: '{:?}'", parameter)
        }

        /// The fasta comment is missing an expected parameter.
        BCalm2MissingParameterError(parameter: String) {
            description("missing parameter")
            display("missing parameter: '{:?}'", parameter)
        }

        /// A node id (from the given graph type) cannot be cast into `usize`.
        BCalm2NodeIdOutOfPrintingRange {
            description("node id is out of range (usize) for displaying")
            display("node id is out of range (usize) for displaying")
        }

        /// A node is not mapped to a mirror node.
        BCalm2NodeWithoutMirror {
            description("node has no mirror")
            display("node has no mirror")
        }

        /// An edge is not mapped to a mirror node.
        BCalm2EdgeWithoutMirror {
            description("edge has no mirror")
            display("edge has no mirror")
        }
    }
}

/// Node data of a bcalm2 node, containing only the data the is typically needed.
#[derive(Debug)]
pub struct BCalm2NodeData {
    // TODO
}

/// The raw node data of a bcalm2 node, including edge information and redundant information (sequence length).
#[derive(Debug, Clone)]
pub struct PlainBCalm2NodeData<GenomeSequenceStoreHandle> {
    /// The numeric id of the bcalm2 node.
    pub id: usize,
    /// The sequence of the bcalm2 node.
    pub sequence_handle: GenomeSequenceStoreHandle,
    /// False if the sequence handle points to the reverse complement of this nodes sequence rather than the actual sequence.
    pub forwards: bool,
    /// The length of the sequence of the bcalm2 node.
    pub length: usize,
    /// The total k-mer abundance of the sequence of the bcalm2 node.
    pub total_abundance: usize,
    /// The mean k-mer abundance of the sequence of the bcalm2 node.
    pub mean_abundance: f64,
    /// The edges stored at the bcalm2 node.
    pub edges: Vec<PlainBCalm2Edge>,
}

/// The raw edge information of a bcalm2 node.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PlainBCalm2Edge {
    /// `true` means `+`, `false` means `-´
    from_side: bool,
    to_node: usize,
    /// `true` means `+`, `false` means `-´
    to_side: bool,
}

impl<GenomeSequenceStoreHandle: Clone> BidirectedData
    for PlainBCalm2NodeData<GenomeSequenceStoreHandle>
{
    fn mirror(&self) -> Self {
        let mut result = self.clone();
        result.forwards = !result.forwards;
        result
    }
}

impl<AlphabetType: Alphabet, GenomeSequenceStore: SequenceStore<AlphabetType>>
    SequenceData<AlphabetType, GenomeSequenceStore>
    for PlainBCalm2NodeData<GenomeSequenceStore::Handle>
{
    fn sequence_handle(&self) -> &GenomeSequenceStore::Handle {
        &self.sequence_handle
    }

    fn sequence_ref<'a>(
        &self,
        source_sequence_store: &'a GenomeSequenceStore,
    ) -> Option<&'a <GenomeSequenceStore as SequenceStore<AlphabetType>>::SequenceRef> {
        if self.forwards {
            let handle = <PlainBCalm2NodeData<GenomeSequenceStore::Handle> as SequenceData<
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
        let handle = <PlainBCalm2NodeData<GenomeSequenceStore::Handle> as SequenceData<
            AlphabetType,
            GenomeSequenceStore,
        >>::sequence_handle(self);

        if self.forwards {
            source_sequence_store.get(handle).convert()
        } else {
            source_sequence_store
                .get(handle)
                .convert_with_reverse_complement()
        }
    }
}

impl<GenomeSequenceStoreHandle: PartialEq> PartialEq
    for PlainBCalm2NodeData<GenomeSequenceStoreHandle>
{
    fn eq(&self, other: &Self) -> bool {
        self.sequence_handle == other.sequence_handle && self.forwards == other.forwards
    }
}

impl<GenomeSequenceStoreHandle: Eq> Eq for PlainBCalm2NodeData<GenomeSequenceStoreHandle> {}

fn parse_bcalm2_fasta_record<
    AlphabetType: Alphabet + 'static,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
>(
    record: Record,
    target_sequence_store: &mut GenomeSequenceStore,
) -> crate::error::Result<PlainBCalm2NodeData<GenomeSequenceStore::Handle>> {
    let id = record
        .id()
        .parse()
        .map_err(|e| Error::with_chain(e, ErrorKind::BCalm2IdError(record.id().to_owned())))?;
    let sequence_handle = target_sequence_store
        .add_from_slice_u8(record.seq())
        .unwrap_or_else(|error| panic!("Genome sequence with id {id} is invalid: {error:?}"));
    let sequence = target_sequence_store.get(&sequence_handle);

    let mut length = None;
    let mut total_abundance = None;
    let mut mean_abundance = None;
    let mut edges = Vec::new();

    for parameter in record.desc().unwrap_or("").split_whitespace() {
        ensure!(
            parameter.len() >= 5,
            Error::from(ErrorKind::BCalm2UnknownParameterError(parameter.to_owned()))
        );
        match &parameter[0..5] {
            "LN:i:" => {
                ensure!(
                    length.is_none(),
                    Error::from(ErrorKind::BCalm2DuplicateParameterError(
                        parameter.to_owned()
                    ))
                );
                length = Some(parameter[5..].parse().map_err(|e| {
                    Error::with_chain(
                        e,
                        ErrorKind::BCalm2MalformedParameterError(parameter.to_owned()),
                    )
                }));
            }
            "KC:i:" => {
                ensure!(
                    total_abundance.is_none(),
                    Error::from(ErrorKind::BCalm2DuplicateParameterError(
                        parameter.to_owned()
                    ))
                );
                total_abundance = Some(parameter[5..].parse().map_err(|e| {
                    Error::with_chain(
                        e,
                        ErrorKind::BCalm2MalformedParameterError(parameter.to_owned()),
                    )
                }));
            }
            "KM:f:" | "km:f:" => {
                ensure!(
                    mean_abundance.is_none(),
                    Error::from(ErrorKind::BCalm2DuplicateParameterError(
                        parameter.to_owned()
                    ))
                );
                mean_abundance = Some(parameter[5..].parse().map_err(|e| {
                    Error::with_chain(
                        e,
                        ErrorKind::BCalm2MalformedParameterError(parameter.to_owned()),
                    )
                }));
            }
            _ => match &parameter[0..2] {
                "L:" => {
                    let parts: Vec<_> = parameter.split(':').collect();
                    ensure!(
                        parts.len() == 4,
                        Error::from(ErrorKind::BCalm2MalformedParameterError(
                            parameter.to_owned()
                        ))
                    );
                    let forward_reverse_to_bool = |c| match c {
                        "+" => Ok(true),
                        "-" => Ok(false),
                        _ => Err(Error::from(ErrorKind::BCalm2MalformedParameterError(
                            parameter.to_owned(),
                        ))),
                    };
                    let from_side = forward_reverse_to_bool(parts[1])?;
                    let to_node = parts[2].parse().map_err(|e| {
                        Error::with_chain(
                            e,
                            ErrorKind::BCalm2MalformedParameterError(parameter.to_owned()),
                        )
                    })?;
                    let to_side = forward_reverse_to_bool(parts[3])?;
                    edges.push(PlainBCalm2Edge {
                        from_side,
                        to_node,
                        to_side,
                    });
                }
                _ => bail!(Error::from(ErrorKind::BCalm2UnknownParameterError(
                    parameter.to_owned()
                ))),
            },
        }
    }

    let length = length.unwrap_or_else(|| {
        bail!(Error::from(ErrorKind::BCalm2MissingParameterError(
            "length (LN)".to_owned()
        )))
    })?;
    ensure!(
        length == sequence.len(),
        Error::from(ErrorKind::BCalm2LengthError(length, sequence.len()))
    );
    let total_abundance = total_abundance.unwrap_or_else(|| {
        bail!(Error::from(ErrorKind::BCalm2MissingParameterError(
            "total abundance (KC)".to_owned()
        )))
    })?;
    let mean_abundance = mean_abundance.unwrap_or_else(|| {
        bail!(Error::from(ErrorKind::BCalm2MissingParameterError(
            "mean abundance (KM)".to_owned()
        )))
    })?;

    Ok(PlainBCalm2NodeData {
        id,
        sequence_handle,
        forwards: true,
        length,
        total_abundance,
        mean_abundance,
        edges,
    })
}

impl<'a, GenomeSequenceStoreHandle: Clone> From<&'a PlainBCalm2NodeData<GenomeSequenceStoreHandle>>
    for PlainBCalm2NodeData<GenomeSequenceStoreHandle>
{
    fn from(data: &'a PlainBCalm2NodeData<GenomeSequenceStoreHandle>) -> Self {
        data.clone()
    }
}

/////////////////////////////
////// NODE CENTRIC IO //////
/////////////////////////////

/// Read a genome graph in bcalm2 fasta format into a node-centric representation from a file.
pub fn read_bigraph_from_bcalm2_as_node_centric_from_file<
    P: AsRef<Path> + Debug,
    AlphabetType: Alphabet + 'static,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
    NodeData: From<PlainBCalm2NodeData<GenomeSequenceStore::Handle>> + BidirectedData,
    EdgeData: Default + Clone,
    Graph: DynamicNodeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    path: P,
    target_sequence_store: &mut GenomeSequenceStore,
) -> crate::error::Result<Graph> {
    read_bigraph_from_bcalm2_as_node_centric(
        bio::io::fasta::Reader::from_file(path).map_err(Error::from)?,
        target_sequence_store,
    )
}

/// Read a genome graph in bcalm2 fasta format into a node-centric representation.
pub fn read_bigraph_from_bcalm2_as_node_centric<
    R: std::io::BufRead,
    AlphabetType: Alphabet + 'static,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
    NodeData: From<PlainBCalm2NodeData<GenomeSequenceStore::Handle>> + BidirectedData,
    EdgeData: Default + Clone,
    Graph: DynamicNodeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    reader: bio::io::fasta::Reader<R>,
    target_sequence_store: &mut GenomeSequenceStore,
) -> crate::error::Result<Graph> {
    struct BiEdge {
        from_node: usize,
        plain_edge: PlainBCalm2Edge,
    }

    let mut bigraph = Graph::default();
    let mut edges = Vec::new();

    for record in reader.records() {
        let record: PlainBCalm2NodeData<GenomeSequenceStore::Handle> =
            parse_bcalm2_fasta_record(record.map_err(Error::from)?, target_sequence_store)?;
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

fn write_plain_bcalm2_node_data_to_bcalm2<GenomeSequenceStoreHandle>(
    node: &PlainBCalm2NodeData<GenomeSequenceStoreHandle>,
    out_neighbors: Vec<(bool, usize, bool)>,
) -> crate::error::Result<String> {
    let mut result = String::new();
    write!(
        result,
        "LN:i:{} KC:i:{} km:f:{:.1}",
        node.length, node.total_abundance, node.mean_abundance
    )
    .map_err(Error::from)?;
    for (node_type, neighbor_id, neighbor_type) in out_neighbors {
        write!(
            result,
            " L:{}:{}:{}",
            if node_type { "+" } else { "-" },
            <usize as NumCast>::from(neighbor_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfPrintingRange))?,
            if neighbor_type { "+" } else { "-" }
        )
        .map_err(Error::from)?;
    }
    Ok(result)
}

/// Write a genome graph in bcalm2 fasta format from a node-centric representation to a file.
pub fn write_node_centric_bigraph_to_bcalm2_to_file<
    P: AsRef<Path>,
    AlphabetType: Alphabet,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
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
    AlphabetType: Alphabet,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
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

/////////////////////////////
////// EDGE CENTRIC IO //////
/////////////////////////////

/// Read a genome graph in bcalm2 fasta format into an edge-centric representation from a file.
pub fn read_bigraph_from_bcalm2_as_edge_centric_from_file<
    P: AsRef<Path> + Debug,
    AlphabetType: Alphabet + 'static + Hash + Eq + Clone,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
    NodeData: Default + Clone,
    EdgeData: From<PlainBCalm2NodeData<GenomeSequenceStore::Handle>> + Clone + Eq + BidirectedData,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    path: P,
    target_sequence_store: &mut GenomeSequenceStore,
    kmer_size: usize,
) -> crate::error::Result<Graph>
where
    <GenomeSequenceStore as SequenceStore<AlphabetType>>::Handle: Clone,
{
    read_bigraph_from_bcalm2_as_edge_centric(
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

/// Read a genome graph in bcalm2 fasta format into an edge-centric representation.
pub fn read_bigraph_from_bcalm2_as_edge_centric<
    R: std::io::BufRead,
    AlphabetType: Alphabet + Hash + Eq + Clone + 'static,
    GenomeSequenceStore: SequenceStore<AlphabetType>,
    NodeData: Default + Clone,
    EdgeData: From<PlainBCalm2NodeData<GenomeSequenceStore::Handle>> + Clone + Eq + BidirectedData,
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
        let record: PlainBCalm2NodeData<GenomeSequenceStore::Handle> =
            parse_bcalm2_fasta_record(record.map_err(Error::from)?, target_sequence_store)?;
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

/// Write a genome graph in bcalm2 fasta format from an edge-centric representation to a file.
pub fn write_edge_centric_bigraph_to_bcalm2_to_file<
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
    for<'a> PlainBCalm2NodeData<GenomeSequenceStore::Handle>: From<&'a EdgeData>,
{
    write_edge_centric_bigraph_to_bcalm2(
        graph,
        source_sequence_store,
        bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

/// Write a genome graph in bcalm2 fasta format from an edge-centric representation.
pub fn write_edge_centric_bigraph_to_bcalm2<
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
    for<'a> PlainBCalm2NodeData<GenomeSequenceStore::Handle>: From<&'a EdgeData>,
{
    let mut output_edges = vec![false; graph.edge_count()];

    for edge_id in graph.edge_indices() {
        if !output_edges[graph
            .mirror_edge_edge_centric(edge_id)
            .ok_or_else(|| Error::from(ErrorKind::BCalm2EdgeWithoutMirror))?
            .as_usize()]
        {
            output_edges[edge_id.as_usize()] = true;
        }
    }

    for edge_id in graph.edge_indices() {
        if output_edges[edge_id.as_usize()] {
            let node_data = PlainBCalm2NodeData::from(graph.edge_data(edge_id));
            let mirror_edge_id = graph
                .mirror_edge_edge_centric(edge_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2EdgeWithoutMirror))?;
            let to_node_plus = graph.edge_endpoints(edge_id).to_node;
            let to_node_minus = graph.edge_endpoints(mirror_edge_id).to_node;

            let mut out_neighbors_plus = Vec::new();
            let mut out_neighbors_minus = Vec::new();

            for neighbor in graph.out_neighbors(to_node_plus) {
                let neighbor_edge_id = neighbor.edge_id.as_usize();

                out_neighbors_plus.push((
                    true,
                    if output_edges[neighbor_edge_id] {
                        PlainBCalm2NodeData::from(graph.edge_data(neighbor.edge_id)).id
                    } else {
                        PlainBCalm2NodeData::from(
                            graph.edge_data(
                                graph
                                    .mirror_edge_edge_centric(neighbor.edge_id)
                                    .ok_or_else(|| {
                                        Error::from(ErrorKind::BCalm2EdgeWithoutMirror)
                                    })?,
                            ),
                        )
                        .id
                    },
                    output_edges[neighbor_edge_id],
                ));
            }
            for neighbor in graph.out_neighbors(to_node_minus) {
                let neighbor_edge_id = neighbor.edge_id.as_usize();

                out_neighbors_minus.push((
                    false,
                    if output_edges[neighbor_edge_id] {
                        PlainBCalm2NodeData::from(graph.edge_data(neighbor.edge_id)).id
                    } else {
                        PlainBCalm2NodeData::from(
                            graph.edge_data(
                                graph
                                    .mirror_edge_edge_centric(neighbor.edge_id)
                                    .ok_or_else(|| {
                                        Error::from(ErrorKind::BCalm2EdgeWithoutMirror)
                                    })?,
                            ),
                        )
                        .id
                    },
                    output_edges[neighbor_edge_id],
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

#[cfg(test)]
mod tests {
    use crate::io::bcalm2::{
        read_bigraph_from_bcalm2_as_edge_centric, read_bigraph_from_bcalm2_as_node_centric,
        write_edge_centric_bigraph_to_bcalm2, write_node_centric_bigraph_to_bcalm2,
    };
    use crate::types::{PetBCalm2EdgeGraph, PetBCalm2NodeGraph};
    use compact_genome::implementation::DefaultSequenceStore;
    use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;

    #[test]
    fn test_node_read_write() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:1:-\n\
            AGT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            GGTCTCGGGTAAGT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:1:-\n\
            ATGATG\n";
        let input = Vec::from(test_file);
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2NodeGraph<_> = read_bigraph_from_bcalm2_as_node_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
        )
        .unwrap();
        let mut output = Vec::new();
        write_node_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:1:-\n\
            AGT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            AATCTCGGGTAAAC\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:1:-\n\
            ACGAGG\n";
        let input = Vec::from(test_file);
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2EdgeGraph<_> = read_bigraph_from_bcalm2_as_edge_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
            3,
        )
        .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_self_loops() {
        let test_file: &'static [u8] =
            b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:0:+ L:+:1:- L:-:0:- L:-:2:+\n\
            AAA\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            GGTCTCGGGTAATT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:0:+ L:-:1:-\n\
            TTGATG\n";
        let input = Vec::from(test_file);
        println!("{}", String::from_utf8(input.clone()).unwrap());
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2EdgeGraph<_> = read_bigraph_from_bcalm2_as_edge_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
            3,
        )
        .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_plus_minus_loop() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:0:- L:+:1:- L:+:2:+\n\
            CAT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:1:- L:+:2:+\n\
            GGTCTCGGGTAAAT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:0:- L:-:1:- L:-:2:+\n\
            ATGATT\n";
        let input = Vec::from(test_file);
        println!("{}", String::from_utf8(input.clone()).unwrap());
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2EdgeGraph<_> = read_bigraph_from_bcalm2_as_edge_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
            3,
        )
        .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_minus_plus_loop() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:1:- L:-:0:+\n\
            ATG\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            GGTCTCGGGTAACA\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:1:-\n\
            CAGATT\n";
        let input = Vec::from(test_file);
        println!("{}", String::from_utf8(input.clone()).unwrap());
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2EdgeGraph<_> = read_bigraph_from_bcalm2_as_edge_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
            3,
        )
        .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_plus_minus_and_minus_plus_loop() {
        let test_file: &'static [u8] =
            b">0 LN:i:4 KC:i:4 km:f:3.0 L:+:0:- L:+:1:- L:+:2:+ L:-:0:+\n\
            CGAT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:1:- L:+:2:+\n\
            GGTCTCGGGTAAAT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:0:- L:-:1:- L:-:2:+\n\
            ATGATG\n";
        let input = Vec::from(test_file);
        println!("{}", String::from_utf8(input.clone()).unwrap());
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2EdgeGraph<_> = read_bigraph_from_bcalm2_as_edge_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
            3,
        )
        .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_all_loops() {
        let test_file: &'static [u8] = b">0 LN:i:4 KC:i:4 km:f:3.0 L:+:0:- L:+:0:+ L:+:1:- L:+:2:+ L:-:0:- L:-:0:+ L:-:1:- L:-:2:+\n\
            ATAT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:0:+ L:+:1:- L:+:2:+\n\
            TCTCGGAAGTAAAT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:0:- L:-:0:+ L:-:1:- L:-:2:+\n\
            ATGATG\n";
        let input = Vec::from(test_file);
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2EdgeGraph<_> = read_bigraph_from_bcalm2_as_edge_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
            3,
        )
        .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_forward_merge() {
        let test_file: &'static [u8] = b"\
            >0 LN:i:3 KC:i:4 km:f:3.0 L:+:2:-\n\
            AGT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:2:-\n\
            GGTCTCGGGTAAGT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:+:0:- L:+:1:-\n\
            AAGAAC\n";
        let input = Vec::from(test_file);
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2EdgeGraph<_> = read_bigraph_from_bcalm2_as_edge_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
            3,
        )
        .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_multigraph() {
        let test_file: &'static [u8] = b"\
            >0 LN:i:7 KC:i:4 km:f:3.0 L:+:2:+ L:-:2:-\n\
            AGTTCTC\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:2:+ L:-:2:-\n\
            AGTCTCGGGTAATC\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:+:0:+ L:+:1:+ L:-:0:- L:-:1:-\n\
            TCGAAG\n";
        let input = Vec::from(test_file);
        let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

        let graph: PetBCalm2EdgeGraph<_> = read_bigraph_from_bcalm2_as_edge_centric(
            bio::io::fasta::Reader::new(test_file),
            &mut sequence_store,
            3,
        )
        .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(
            &graph,
            &sequence_store,
            bio::io::fasta::Writer::new(&mut output),
        )
        .unwrap();

        debug_assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }
}
