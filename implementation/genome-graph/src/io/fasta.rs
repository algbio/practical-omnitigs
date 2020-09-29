use crate::io::bcalm2::PlainBCalm2NodeData;
use bigraph::traitgraph::interface::ImmutableGraphContainer;
use bigraph::traitgraph::walks::EdgeWalk;
use compact_genome::interface::{ExtendableGenome, Genome};
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

/// Write a sequence of walks in a graph as fasta records.
pub fn write_walks_as_fasta<
    'ws,
    Graph: ImmutableGraphContainer<EdgeData = PlainBCalm2NodeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph>,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    Writer: std::io::Write,
>(
    graph: &Graph,
    kmer_size: usize,
    walks: WalkSource,
    writer: &mut bio::io::fasta::Writer<Writer>,
) -> crate::error::Result<()> {
    for (i, walk) in walks.into_iter().enumerate() {
        if walk.is_empty() {
            return Err(Error::from_kind(ErrorKind::EmptyWalkError).into());
        }

        let mut sequence = graph.edge_data(walk[0]).sequence.clone();
        for edge in walk.iter().skip(1) {
            let edge_sequence = &graph.edge_data(*edge).sequence[kmer_size - 1..];
            sequence.extend(edge_sequence.iter().copied());
        }

        let record =
            bio::io::fasta::Record::with_attrs(&format!("{}", i), None, &sequence.into_vec());
        writer.write_record(&record).map_err(Error::from)?;
    }

    Ok(())
}

/// Write a sequence of walks in a graph as fasta records to a file.
/// The given file is created if it does not exist or truncated if it does exist.
pub fn write_walks_as_fasta_file<
    'ws,
    Graph: ImmutableGraphContainer<EdgeData = PlainBCalm2NodeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph>,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
    P: AsRef<Path>,
>(
    graph: &Graph,
    kmer_size: usize,
    walks: WalkSource,
    path: P,
) -> crate::error::Result<()> {
    write_walks_as_fasta(
        graph,
        kmer_size,
        walks,
        &mut bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}
