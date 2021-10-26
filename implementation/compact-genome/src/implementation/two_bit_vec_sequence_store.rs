//! A vector sequence store where each character is encoded as two bits.

use crate::implementation::two_bit_vec_sequence::{TwoBitVectorGenome, TwoBitVectorSubGenome};
use crate::interface::sequence::GenomeSequence;
use crate::interface::sequence_store::{
    HandleWithLength, HandleWithSubsequence, InverseMappingSequenceStore, SequenceStore,
};
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

impl InverseMappingSequenceStore for TwoBitVectorSequenceStore {
    fn map_sequence_ref_to_handle(&self, sequence_ref: &Self::SequenceRef) -> Self::Handle {
        let raw_offset = self.sequence.bits.offset_from(&sequence_ref.bits[..]);
        debug_assert!(raw_offset >= 0);
        let offset = raw_offset as usize / 2;

        Self::Handle {
            offset,
            len: sequence_ref.len(),
        }
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

#[cfg(test)]
mod tests {
    use crate::implementation::two_bit_vec_sequence_store::TwoBitVectorSequenceStore;
    use crate::interface::sequence::GenomeSequence;
    use crate::interface::sequence_store::{InverseMappingSequenceStore, SequenceStore};

    #[test]
    fn test_inverse_mapping() {
        let mut sequence_store = TwoBitVectorSequenceStore::new();
        let handle1 = sequence_store.add(&b"ACGTTG"[..]);
        let handle2 = sequence_store.add(&b"CGACTG"[..]);
        let reference1 = sequence_store.get(&handle1);
        let reference2 = sequence_store.get(&handle2);
        debug_assert_eq!(&reference1.convert::<Vec<_>, _>(), b"ACGTTG");
        debug_assert_eq!(&reference2.convert::<Vec<_>, _>(), b"CGACTG");
        debug_assert_eq!(
            sequence_store.map_sequence_ref_to_handle(reference1),
            handle1
        );
        debug_assert_eq!(
            sequence_store.map_sequence_ref_to_handle(reference2),
            handle2
        );
    }
}
