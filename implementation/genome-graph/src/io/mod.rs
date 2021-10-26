use compact_genome::implementation::ascii_vec_sequence_store::{AsciiVectorSequenceStore, AsciiVectorSequenceStoreHandle};
use compact_genome::implementation::two_bit_vec_sequence_store::{TwoBitVectorSequenceStore, TwoBitVectorSequenceStoreHandle};
use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
use compact_genome::interface::sequence_store::SequenceStore;

/// A module providing types and functions for IO in the bcalm2 fasta format.
pub mod bcalm2;
/// A module providing functions to read and write walks in a de Bruijn graph as fasta.
pub mod fasta;
/// A module providing types and functions for IO in gfa format.
pub mod gfa;
/// A module providing types and functions for IO in the wtdbg2 graph and contig formats.
pub mod wtdbg2;

/// Node or edge data of a genome graph that has an associated sequence.
pub trait SequenceData<GenomeSequenceStore: SequenceStore> {
    /// Returns the handle of the sequence stored in this type.
    fn sequence_handle(&self) -> &GenomeSequenceStore::Handle;

    /// Returns a reference to the sequence pointed to by the handle of this type.
    /// If this type has a sequence that is not stored anywhere (e.g. a reverse complement), it returns `None`.
    fn sequence_ref<'a>(
        &self,
        source_sequence_store: &'a GenomeSequenceStore,
    ) -> Option<&'a GenomeSequenceStore::SequenceRef>;

    /// Returns an owned copy of the sequence pointed to by the handle of this type.
    fn sequence_owned<
        ResultSequence: for<'a> OwnedGenomeSequence<'a, ResultSubsequence>,
        ResultSubsequence: for<'a> GenomeSequence<'a, ResultSubsequence> + ?Sized,
    >(
        &self,
        source_sequence_store: &GenomeSequenceStore,
    ) -> ResultSequence;
}

impl SequenceData<TwoBitVectorSequenceStore> for TwoBitVectorSequenceStoreHandle {
    fn sequence_handle(&self) -> &<TwoBitVectorSequenceStore as SequenceStore>::Handle {
        &self
    }

    fn sequence_ref<'a>(&self, source_sequence_store: &'a TwoBitVectorSequenceStore) -> Option<&'a <TwoBitVectorSequenceStore as SequenceStore>::SequenceRef> {
        Some(source_sequence_store.get(self))
    }

    fn sequence_owned<ResultSequence: for<'a> OwnedGenomeSequence<'a, ResultSubsequence>, ResultSubsequence: for<'a> GenomeSequence<'a, ResultSubsequence> + ?Sized>(&self, source_sequence_store: &TwoBitVectorSequenceStore) -> ResultSequence {
        source_sequence_store.get(self).convert()
    }
}

impl SequenceData<AsciiVectorSequenceStore> for AsciiVectorSequenceStoreHandle {
    fn sequence_handle(&self) -> &<AsciiVectorSequenceStore as SequenceStore>::Handle {
        &self
    }

    fn sequence_ref<'a>(&self, source_sequence_store: &'a AsciiVectorSequenceStore) -> Option<&'a <AsciiVectorSequenceStore as SequenceStore>::SequenceRef> {
        Some(source_sequence_store.get(self))
    }

    fn sequence_owned<ResultSequence: for<'a> OwnedGenomeSequence<'a, ResultSubsequence>, ResultSubsequence: for<'a> GenomeSequence<'a, ResultSubsequence> + ?Sized>(&self, source_sequence_store: &AsciiVectorSequenceStore) -> ResultSequence {
        source_sequence_store.get(self).convert()
    }
}