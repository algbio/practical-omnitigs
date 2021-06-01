use crate::io::SequenceData;
use bigraph::traitgraph::interface::ImmutableGraphContainer;
use bigraph::traitgraph::traitsequence::interface::Sequence;
use bigraph::traitgraph::walks::{EdgeWalk, NodeWalk};
use compact_genome::implementation::DefaultGenome;
use compact_genome::interface::sequence::{
    EditableGenomeSequence, GenomeSequence, OwnedGenomeSequence,
};
use compact_genome::interface::sequence_store::SequenceStore;
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
    }

    errors {
        /// A walk is empty.
        EmptyWalkError {
            description("walk is empty")
            display("walk is empty")
        }
    }
}

/// Data that can be output as fasta record.
pub trait FastaData<SourceSequenceStore: SequenceStore> {
    /// The type storing the genome sequence of this fasta record.
    type Genome: for<'a> OwnedGenomeSequence<'a, Self::GenomeSubsequence>
        + for<'a> EditableGenomeSequence<'a, Self::GenomeSubsequence>;
    /// The subsequence type of `Genome`.
    type GenomeSubsequence: for<'a> GenomeSequence<'a, Self::GenomeSubsequence> + ?Sized;

    /// Returns the sequence of this fasta record.
    fn sequence<'a>(
        &self,
        source_sequence_store: &'a SourceSequenceStore,
    ) -> &'a Self::GenomeSubsequence;
}

/// Write a sequence of walks in a graph as fasta records.
pub fn write_walks_as_fasta<
    'ws,
    SourceSequenceStore: SequenceStore,
    EdgeData: SequenceData<SourceSequenceStore>,
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

        let mut sequence: DefaultGenome = graph
            .edge_data(walk[0])
            .sequence_owned(source_sequence_store);
        for edge in walk.iter().skip(1) {
            let edge_data = graph.edge_data(*edge);
            if let Some(sequence_ref) = edge_data.sequence_ref(source_sequence_store) {
                let sequence_ref = sequence_ref.iter().skip(kmer_size - 1);
                sequence.extend(sequence_ref.copied());
            } else {
                let sequence_owned: DefaultGenome = edge_data.sequence_owned(source_sequence_store);
                let sequence_owned = sequence_owned.iter().skip(kmer_size - 1);
                sequence.extend(sequence_owned.copied());
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
    SourceSequenceStore: SequenceStore,
    EdgeData: SequenceData<SourceSequenceStore>,
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
    SourceSequenceStore: SequenceStore,
    NodeData: SequenceData<SourceSequenceStore>,
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

        let mut sequence: DefaultGenome = graph
            .node_data(walk[0])
            .sequence_owned(source_sequence_store);
        for node in walk.iter().skip(1) {
            let node_data = graph.node_data(*node);
            if let Some(sequence_ref) = node_data.sequence_ref(source_sequence_store) {
                let sequence_ref = sequence_ref.iter().skip(kmer_size - 1);
                sequence.extend(sequence_ref.copied());
            } else {
                let sequence_owned: DefaultGenome = node_data.sequence_owned(source_sequence_store);
                let sequence_owned = sequence_owned.iter().skip(kmer_size - 1);
                sequence.extend(sequence_owned.copied());
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
    SourceSequenceStore: SequenceStore,
    NodeData: SequenceData<SourceSequenceStore>,
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
