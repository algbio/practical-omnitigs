//! A vector sequence store where each character is encoded as two bits.

use crate::implementation::two_bit_vec_sequence::{TwoBitVectorGenome, TwoBitVectorSubGenome};
use crate::interface::sequence::GenomeSequence;
use crate::interface::sequence_store::{HandleWithLength, SequenceStore};
use traitsequence::interface::Sequence;

/// A bitvector based sequence store.
#[derive(Default, Clone, Eq, PartialEq, Debug)]
pub struct TwoBitVectorSequenceStore {
    sequence: TwoBitVectorGenome,
}

/// A handle of a sequence in an [TwoBitVectorSequenceStore].
#[derive(Default, Debug, Clone, Eq, PartialEq)]
pub struct TwoBitVectorSequenceStoreHandle {
    offset: usize,
    len: usize,
}

impl TwoBitVectorSequenceStore {
    /// Creates a new instance.
    pub fn new() -> Self {
        Self {
            sequence: Default::default(),
        }
    }
}

impl SequenceStore for TwoBitVectorSequenceStore {
    type Handle = TwoBitVectorSequenceStoreHandle;
    type SequenceRef = TwoBitVectorSubGenome;

    fn add<
        Sequence: for<'a> GenomeSequence<'a, Subsequence> + ?Sized,
        Subsequence: for<'a> GenomeSequence<'a, Subsequence> + ?Sized,
    >(
        &mut self,
        s: &Sequence,
    ) -> Self::Handle {
        let offset = self.sequence.len();
        let len = s.len();
        self.sequence.extend(s.iter().copied());
        Self::Handle { offset, len }
    }

    fn get(&self, handle: &Self::Handle) -> &Self::SequenceRef {
        &self.sequence[handle.offset..handle.offset + handle.len]
    }
}

impl HandleWithLength for TwoBitVectorSequenceStoreHandle {
    fn len(&self) -> usize {
        self.len
    }
}
