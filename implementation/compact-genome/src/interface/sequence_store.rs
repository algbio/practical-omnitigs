//! Traits for genome sequence stores.

use crate::interface::sequence::{GenomeSequence, OwnedGenomeSequence};

/// A store for sequence data, which gives out handles that can be used to retrieve concrete sequences.
pub trait SequenceStore {
    /// A handle to a sequence in this store. Can be used to retrieve the respective sequence.
    type Handle;

    /// A reference to a sequence stored in this store.
    /// This trait only uses `&SequenceRef`, so `SequenceRef` should be the base type and not a reference type itself.
    type SequenceRef: for<'a> GenomeSequence<'a, Self::SequenceRef> + ?Sized;

    /// Adds a sequence to this store and returns a handle for later retrieval.
    /// Handles do not borrow the sequence store, so they can exist while the store is modified.
    fn add<
        Sequence: for<'a> OwnedGenomeSequence<'a, Subsequence>,
        Subsequence: for<'a> GenomeSequence<'a, Subsequence>,
    >(
        &mut self,
        s: Sequence,
    ) -> Self::Handle;

    /// Returns a reference to a sequence in this store, given the handle.
    /// The reference borrows the sequence store, so it cannot be mutated while references exist.
    /// On the other hand, handles do not borrow, so they can exist while the store is modified.
    fn get(&self, handle: &Self::Handle) -> &Self::SequenceRef;
}
