use bigraph::traitgraph::interface::ImmutableGraphContainer;
use bigraph::traitgraph::traitsequence::interface::Sequence;
use bigraph::traitgraph::walks::{EdgeWalk, NodeWalk};
use std::path::Path;
use compact_genome::interface::{GenomeSequence, EditableGenomeSequence, OwnedGenomeSequence};
use compact_genome::implementation::vec::AsciiVectorGenome;
use std::iter::FromIterator;

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
pub trait FastaData
{
    /// The type storing the genome sequence of this fasta record.
    type Genome: for<'a> OwnedGenomeSequence<'a, Self::GenomeSubsequence> + for<'a> EditableGenomeSequence<'a, Self::GenomeSubsequence>;
    /// The subsequence type of `Genome`.
    type GenomeSubsequence: for <'a> GenomeSequence<'a, Self::GenomeSubsequence> + ?Sized;

    /// Returns the sequence of this fasta record.
    fn sequence(&self) -> &Self::Genome;
}

/// Write a sequence of walks in a graph as fasta records.
pub fn write_walks_as_fasta<
    'ws,
    EdgeData: FastaData,
    Graph: ImmutableGraphContainer<EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    Writer: std::io::Write,
>(
    graph: &Graph,
    kmer_size: usize,
    walks: WalkSource,
    writer: &mut bio::io::fasta::Writer<Writer>,
) -> crate::error::Result<()>
{
    for (i, walk) in walks.into_iter().enumerate() {
        if walk.is_empty() {
            return Err(Error::from_kind(ErrorKind::EmptyWalkError).into());
        }

        let mut sequence = AsciiVectorGenome::from_iter(graph.edge_data(walk[0]).sequence().iter().copied());
        for edge in walk.iter().skip(1) {
            let edge_sequence = graph.edge_data(*edge).sequence();
            let edge_sequence = edge_sequence.iter().skip(kmer_size - 1);
            sequence.extend(edge_sequence.copied());
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
    EdgeData: FastaData,
    Graph: ImmutableGraphContainer<EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    P: AsRef<Path>,
>(
    graph: &Graph,
    kmer_size: usize,
    walks: WalkSource,
    path: P,
) -> crate::error::Result<()>
{
    write_walks_as_fasta(
        graph,
        kmer_size,
        walks,
        &mut bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

/// Write a sequence of node-centric walks in a graph as fasta records.
pub fn write_node_centric_walks_as_fasta<
    'ws,
    NodeData: FastaData,
    Graph: ImmutableGraphContainer<NodeData = NodeData>,
    Walk: 'ws + for<'w> NodeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> NodeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    Writer: std::io::Write,
>(
    graph: &Graph,
    kmer_size: usize,
    walks: WalkSource,
    writer: &mut bio::io::fasta::Writer<Writer>,
) -> crate::error::Result<()>
{
    for (i, walk) in walks.into_iter().enumerate() {
        if walk.is_empty() {
            return Err(Error::from_kind(ErrorKind::EmptyWalkError).into());
        }

        let mut sequence = AsciiVectorGenome::from_iter(graph.node_data(walk[0]).sequence().iter().copied());
        for node in walk.iter().skip(1) {
            let node_sequence = graph.node_data(*node).sequence();
            let node_sequence = node_sequence.iter().skip(kmer_size - 1);
            sequence.extend(node_sequence.copied());
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
    NodeData: FastaData,
    Graph: ImmutableGraphContainer<NodeData = NodeData>,
    Walk: 'ws + for<'w> NodeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> NodeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    P: AsRef<Path>,
>(
    graph: &Graph,
    kmer_size: usize,
    walks: WalkSource,
    path: P,
) -> crate::error::Result<()>
{
    write_node_centric_walks_as_fasta(
        graph,
        kmer_size,
        walks,
        &mut bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}
