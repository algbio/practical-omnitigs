//! An ASCII vector based sequence store.

use crate::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
use crate::interface::sequence_store::SequenceStore;

/// An ASCII vector based sequence store.
pub struct AsciiVectorSequenceStore {
    sequence: Vec<u8>,
}

/// A handle of a sequence in an [AsciiVectorSequenceStore].
pub struct AsciiVectorSequenceStoreHandle {
    offset: usize,
    len: usize,
}

impl AsciiVectorSequenceStore {
    /// Creates a new instance.
    pub fn new() -> Self {
        Self {
            sequence: Vec::new(),
        }
    }
}

impl SequenceStore for AsciiVectorSequenceStore {
    type Handle = AsciiVectorSequenceStoreHandle;
    type SequenceRef = [u8];

    fn add<
        Sequence: for<'a> OwnedGenomeSequence<'a, Subsequence>,
        Subsequence: for<'a> GenomeSequence<'a, Subsequence>,
    >(
        &mut self,
        s: Sequence,
    ) -> Self::Handle {
        let offset = self.sequence.len();
        let len = s.len();
        self.sequence.extend(s.iter());
        Self::Handle { offset, len }
    }

    fn get(&self, handle: &Self::Handle) -> &Self::SequenceRef {
        &self.sequence[handle.offset..handle.offset + handle.len]
    }
}
