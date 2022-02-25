//! Traits for genome sequence stores.

use crate::interface::alphabet::{Alphabet, AlphabetError};
use crate::interface::sequence::GenomeSequence;

/// A store for sequence data, which gives out handles that can be used to retrieve concrete sequences.
pub trait SequenceStore<AlphabetType: Alphabet> {
    /// A handle to a sequence in this store. Can be used to retrieve the respective sequence.
    type Handle;

    /// A reference to a sequence stored in this store.
    /// This trait only uses `&SequenceRef`, so `SequenceRef` should be the base type and not a reference type itself.
    type SequenceRef: for<'a> GenomeSequence<'a, AlphabetType, Self::SequenceRef> + ?Sized;

    /// Adds a sequence to this store and returns a handle for later retrieval.
    /// Handles do not borrow the sequence store, so they can exist while the store is modified.
    fn add<
        Sequence: for<'a> GenomeSequence<'a, AlphabetType, Subsequence> + ?Sized,
        Subsequence: for<'a> GenomeSequence<'a, AlphabetType, Subsequence> + ?Sized,
    >(
        &mut self,
        s: &Sequence,
    ) -> Self::Handle;

    /// Adds a sequence to this store and returns a handle for later retrieval.
    /// Handles do not borrow the sequence store, so they can exist while the store is modified.
    ///
    /// This method expects an `IntoIter` over ASCII characters, and returns `None` if any of the characters is not part of the alphabet.
    fn add_from_iter_u8<IteratorType: IntoIterator<Item = u8>>(
        &mut self,
        iter: IteratorType,
    ) -> Result<Self::Handle, AlphabetError>;

    /// Adds a sequence to this store and returns a handle for later retrieval.
    /// Handles do not borrow the sequence store, so they can exist while the store is modified.
    ///
    /// This method expects slice of ASCII characters, and returns `None` if any of the characters is not part of the alphabet.
    fn add_from_slice_u8(&mut self, slice: &[u8]) -> Result<Self::Handle, AlphabetError> {
        self.add_from_iter_u8(slice.iter().copied())
    }

    /// Returns a reference to a sequence in this store, given the handle.
    /// The reference borrows the sequence store, so it cannot be mutated while references exist.
    /// On the other hand, handles do not borrow, so they can exist while the store is modified.
    fn get(&self, handle: &Self::Handle) -> &Self::SequenceRef;
}

/// A sequence store that is able to map from references to sequences back to handles.
pub trait InverseMappingSequenceStore<AlphabetType: Alphabet>: SequenceStore<AlphabetType> {
    /// Returns a handle that refers the given sequence reference.
    fn map_sequence_ref_to_handle(&self, sequence_ref: &Self::SequenceRef) -> Self::Handle;
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
