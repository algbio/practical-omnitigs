//! Traits for genome sequences.

use crate::interface::alphabet::{Alphabet, AlphabetCharacter, AlphabetError};
use std::cmp::Ordering;
use std::iter::{FromIterator, Map, Rev};
use traitsequence::interface::{EditableSequence, OwnedSequence, Sequence, SequenceMut};

pub mod neighbor_iterators;

/// An iterator over the reverse complement of a genome sequence.
pub type ReverseComplementIterator<I, AlphabetType> = Map<
    Rev<I>,
    for<'c> fn(
        &'c <AlphabetType as Alphabet>::CharacterType,
    ) -> <AlphabetType as Alphabet>::CharacterType,
>;

/// A genome sequence.
pub trait GenomeSequence<
    'a,
    AlphabetType: Alphabet,
    GenomeSubsequence: GenomeSequence<'a, AlphabetType, GenomeSubsequence> + ?Sized,
>: Sequence<'a, AlphabetType::CharacterType, GenomeSubsequence>
{
    /// Returns true if this genome is valid, i.e. it contains no invalid characters.
    fn is_valid(&'a self) -> bool {
        true
    }

    /// Copies this genome string into a `Vec`.
    fn clone_as_vec(&'a self) -> Vec<u8> {
        self.iter()
            .cloned()
            .map(AlphabetType::character_to_ascii)
            .collect()
    }

    /// Get a reference to this genome as its subsequence type.
    fn as_genome_subsequence(&self) -> &GenomeSubsequence {
        self.index(0..self.len())
    }

    /// Returns the genome as nucleotide string.
    fn as_string(&'a self) -> String {
        String::from_utf8(self.clone_as_vec())
            .expect("Genome contains non-utf8 characters (It should be ASCII only).")
    }

    /// Returns an iterator over the reverse complement of this genome.
    /// Panics if the iterator his an invalid character (see [not valid](GenomeSequence::is_valid)).
    fn reverse_complement_iter(
        &'a self,
    ) -> ReverseComplementIterator<Self::Iterator, AlphabetType> {
        self.iter()
            .rev()
            .map(AlphabetType::CharacterType::complement)
    }

    /// Returns an owned copy of the reverse complement of this genome.
    /// Panics if this genome is [not valid](GenomeSequence::is_valid).
    fn convert_with_reverse_complement<
        ReverseComplementSequence: for<'rc> OwnedGenomeSequence<'rc, AlphabetType, ReverseComplementSubsequence>,
        ReverseComplementSubsequence: for<'rc> GenomeSequence<'a, AlphabetType, ReverseComplementSubsequence> + ?Sized,
    >(
        &'a self,
    ) -> ReverseComplementSequence {
        self.reverse_complement_iter().collect()
    }

    /// Returns an owned copy of this genome.
    fn convert<
        ResultSequence: for<'rc> OwnedGenomeSequence<'rc, AlphabetType, ResultSubsequence>,
        ResultSubsequence: for<'rc> GenomeSequence<'rc, AlphabetType, ResultSubsequence> + ?Sized,
    >(
        &'a self,
    ) -> ResultSequence {
        self.iter().cloned().collect()
    }

    /// Returns true if the genome is canonical.
    /// A canonical genome is lexicographically smaller or equal to its reverse complement.
    fn is_canonical(&'a self) -> bool {
        for (forward_character, reverse_character) in
            self.iter().cloned().zip(self.reverse_complement_iter())
        {
            match forward_character.cmp(&reverse_character) {
                Ordering::Less => return true,
                Ordering::Greater => return false,
                _ => {}
            }
        }
        true
    }

    /// Returns true if the genome is self-complemental.
    /// A self-complemental genome is equivalent to its reverse complement.
    fn is_self_complemental(&'a self) -> bool {
        self.iter().cloned().eq(self.reverse_complement_iter())
    }
}

/// A genome sequence that is owned, i.e. not a reference.
pub trait OwnedGenomeSequence<
    'a,
    AlphabetType: Alphabet,
    GenomeSubsequence: GenomeSequence<'a, AlphabetType, GenomeSubsequence> + ?Sized,
>:
    for<'s> GenomeSequence<'s, AlphabetType, GenomeSubsequence>
    + FromIterator<AlphabetType::CharacterType>
    + for<'s> OwnedSequence<'s, AlphabetType::CharacterType, GenomeSubsequence>
{
    /// Returns the reverse complement of this genome.
    /// Panics if this genome is [not valid](GenomeSequence::is_valid).
    fn clone_as_reverse_complement(&'a self) -> Self {
        self.reverse_complement_iter().collect()
    }

    /// Constructs an owned genome sequence from an `IntoIter` over ASCII characters.
    /// If any character is not part of the alphabet, then `None` is returned.
    fn from_iter_u8<T: IntoIterator<Item = u8>>(iter: T) -> Result<Self, AlphabetError> {
        iter.into_iter()
            .map(AlphabetType::ascii_to_character)
            .collect()
    }

    /// Constructs an owned genome sequence from a slice of ASCII characters.
    /// If any character is not part of the alphabet, then `None` is returned.
    fn from_slice_u8(slice: &[u8]) -> Result<Self, AlphabetError> {
        Self::from_iter_u8(slice.iter().copied())
    }
}

/// A mutable genome sequence.
pub trait GenomeSequenceMut<
    'a,
    AlphabetType: Alphabet,
    GenomeSubsequenceMut: GenomeSequenceMut<'a, AlphabetType, GenomeSubsequenceMut> + ?Sized,
>:
    SequenceMut<'a, AlphabetType::CharacterType, GenomeSubsequenceMut>
    + GenomeSequence<'a, AlphabetType, GenomeSubsequenceMut>
{
}

type IntoIterU8<SourceType, AlphabetType> = Map<
    <SourceType as IntoIterator>::IntoIter,
    fn(<AlphabetType as Alphabet>::CharacterType) -> u8,
>;

/// An editable genome sequence.
pub trait EditableGenomeSequence<
    'a,
    AlphabetType: Alphabet,
    GenomeSubsequence: GenomeSequence<'a, AlphabetType, GenomeSubsequence> + ?Sized,
>:
    EditableSequence<'a, AlphabetType::CharacterType, GenomeSubsequence>
    + GenomeSequence<'a, AlphabetType, GenomeSubsequence>
{
    /// Converts this genome sequence into an iterator over ASCII characters.
    fn into_iter_u8(self) -> IntoIterU8<Self, AlphabetType> {
        self.into_iter().map(AlphabetType::character_to_ascii)
    }

    /// Extends this genome from a sequence of ASCII characters.
    fn extend_from_iter_u8<IteratorType: IntoIterator<Item = u8>>(
        &mut self,
        iter: IteratorType,
    ) -> Result<(), AlphabetError> {
        let original_len = self.len();
        let iter = iter.into_iter();
        let (size, _) = iter.size_hint();
        self.reserve(size);
        for item in iter {
            match AlphabetType::ascii_to_character(item) {
                Ok(character) => self.push(character),
                Err(error) => {
                    self.resize(
                        original_len,
                        AlphabetType::CharacterType::from_index(0).unwrap(),
                    );
                    return Err(error);
                }
            }
        }

        Ok(())
    }

    /// Extends this genome from a sequence of ASCII characters.
    fn extend_from_slice_u8(&mut self, slice: &[u8]) -> Result<(), AlphabetError> {
        self.extend_from_iter_u8(slice.iter().copied())
    }

    /// Reserve memory for at least `additional` items.
    fn reserve(&mut self, additional: usize);

    /// Resize to contain the given number of items.
    /// Empty spaces are filled with the given default item.
    fn resize(&mut self, len: usize, default: AlphabetType::CharacterType);

    /// Insert the given character at the end of the genome sequence.
    fn push(&mut self, character: AlphabetType::CharacterType);
}
