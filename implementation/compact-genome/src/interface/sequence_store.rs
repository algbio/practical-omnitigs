//! Traits for genome sequence stores.

use crate::interface::sequence::GenomeSequence;

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
        Sequence: for<'a> GenomeSequence<'a, Subsequence> + ?Sized,
        Subsequence: for<'a> GenomeSequence<'a, Subsequence> + ?Sized,
    >(
        &mut self,
        s: &Sequence,
    ) -> Self::Handle;

    /// Returns a reference to a sequence in this store, given the handle.
    /// The reference borrows the sequence store, so it cannot be mutated while references exist.
    /// On the other hand, handles do not borrow, so they can exist while the store is modified.
    fn get(&self, handle: &Self::Handle) -> &Self::SequenceRef;
}

/// A handle of a sequence store that can compute the length of the referred sequence without retrieving the sequence.
pub trait HandleWithLength {
    /// Returns the length of the sequence referred by this handle without retrieving the sequence.
    fn len(&self) -> usize;

    /// Returns true if the length of the sequence referred by this handle is zero, without retrieving the sequence.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A handle that allows the creation of handles referring arbitrary subsequences of this handle.
pub trait HandleWithSubsequence<RangeType> {
    /// Returns a new handle that refers the subsequence of this handle as indicated by the range.
    /// This method may panic if the range does not define a subsequence of the sequence referred by this handle.
    /// However, since the handle might not know anything about the sequence it refers, it might also silently ignore such errors.
    fn subsequence_handle(&self, range: RangeType) -> Self;
}
