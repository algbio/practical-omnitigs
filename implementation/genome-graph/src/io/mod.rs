use compact_genome::implementation::bit_vec_sequence_store::{
    BitVectorSequenceStore, BitVectorSequenceStoreHandle,
};
use compact_genome::implementation::vec_sequence_store::{
    VectorSequenceStore, VectorSequenceStoreHandle,
};
use compact_genome::interface::alphabet::Alphabet;
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
pub trait SequenceData<AlphabetType: Alphabet, GenomeSequenceStore: SequenceStore<AlphabetType>> {
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
        ResultSequence: for<'a> OwnedGenomeSequence<'a, AlphabetType, ResultSubsequence>,
        ResultSubsequence: for<'a> GenomeSequence<'a, AlphabetType, ResultSubsequence> + ?Sized,
    >(
        &self,
        source_sequence_store: &GenomeSequenceStore,
    ) -> ResultSequence;
}

impl<AlphabetType: Alphabet + 'static>
    SequenceData<AlphabetType, BitVectorSequenceStore<AlphabetType>>
    for BitVectorSequenceStoreHandle<AlphabetType>
{
    fn sequence_handle(
        &self,
    ) -> &<BitVectorSequenceStore<AlphabetType> as SequenceStore<AlphabetType>>::Handle {
        self
    }

    fn sequence_ref<'a>(
        &self,
        source_sequence_store: &'a BitVectorSequenceStore<AlphabetType>,
    ) -> Option<
        &'a <BitVectorSequenceStore<AlphabetType> as SequenceStore<AlphabetType>>::SequenceRef,
    > {
        Some(source_sequence_store.get(self))
    }

    fn sequence_owned<
        ResultSequence: for<'a> OwnedGenomeSequence<'a, AlphabetType, ResultSubsequence>,
        ResultSubsequence: for<'a> GenomeSequence<'a, AlphabetType, ResultSubsequence> + ?Sized,
    >(
        &self,
        source_sequence_store: &BitVectorSequenceStore<AlphabetType>,
    ) -> ResultSequence {
        source_sequence_store.get(self).convert()
    }
}

impl<AlphabetType: Alphabet + 'static> SequenceData<AlphabetType, VectorSequenceStore<AlphabetType>>
    for VectorSequenceStoreHandle<AlphabetType>
{
    fn sequence_handle(
        &self,
    ) -> &<VectorSequenceStore<AlphabetType> as SequenceStore<AlphabetType>>::Handle {
        self
    }

    fn sequence_ref<'a>(
        &self,
        source_sequence_store: &'a VectorSequenceStore<AlphabetType>,
    ) -> Option<&'a <VectorSequenceStore<AlphabetType> as SequenceStore<AlphabetType>>::SequenceRef>
    {
        Some(source_sequence_store.get(self))
    }

    fn sequence_owned<
        ResultSequence: for<'a> OwnedGenomeSequence<'a, AlphabetType, ResultSubsequence>,
        ResultSubsequence: for<'a> GenomeSequence<'a, AlphabetType, ResultSubsequence> + ?Sized,
    >(
        &self,
        source_sequence_store: &VectorSequenceStore<AlphabetType>,
    ) -> ResultSequence {
        source_sequence_store.get(self).convert()
    }
}
