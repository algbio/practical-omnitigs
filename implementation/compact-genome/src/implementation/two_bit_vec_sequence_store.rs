//! A vector sequence store where each character is encoded as two bits.

use crate::implementation::two_bit_vec_sequence::{TwoBitVectorGenome, TwoBitVectorSubGenome};
use crate::interface::sequence::GenomeSequence;
use crate::interface::sequence_store::{HandleWithLength, HandleWithSubsequence, SequenceStore};
use traitsequence::interface::Sequence;

/// A bitvector based sequence store.
#[derive(Default, Clone, Eq, PartialEq, Debug)]
pub struct TwoBitVectorSequenceStore {
    sequence: TwoBitVectorGenome,
}

/// A handle of a sequence in an [TwoBitVectorSequenceStore].
#[derive(Default, Debug, Clone, Copy, Eq, PartialEq)]
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

    /// Returns the number of bytes consumed by the characters stored in this sequence store.
    pub fn size_in_memory(&self) -> usize {
        (self.sequence.len() - 1) / 4 + 1 // Rounding up integer division
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

impl HandleWithSubsequence<core::ops::Range<usize>> for TwoBitVectorSequenceStoreHandle {
    fn subsequence_handle(&self, range: core::ops::Range<usize>) -> Self {
        let result = Self {
            offset: self.offset + range.start,
            len: range.end - range.start,
        };
        debug_assert!(self.offset + self.len >= result.offset + result.len);
        result
    }
}
