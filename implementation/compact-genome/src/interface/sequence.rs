//! Traits for genome sequences.

use itertools::Itertools;
use std::iter::{Copied, FromIterator, Map, Rev};
use traitsequence::interface::{EditableSequence, Sequence, SequenceMut};

/// An iterator over the reverse complement of a genome sequence.
pub type ReverseComplementIterator<I> =
    Map<Map<Rev<Copied<I>>, fn(u8) -> Option<u8>>, fn(Option<u8>) -> u8>;

/// A genome sequence.
pub trait GenomeSequence<'a, GenomeSubsequence: GenomeSequence<'a, GenomeSubsequence> + ?Sized>:
    Sequence<'a, u8, GenomeSubsequence>
{
    /// Returns true if this genome is valid, i.e. it contains no invalid characters.
    /// Valid characters are defined by [is_valid_ascii_genome_character()](is_valid_ascii_genome_character)
    fn is_valid(&'a self) -> bool {
        self.iter().copied().all(is_valid_ascii_genome_character)
    }

    /// Returns a duplicate-free vector of all invalid characters in this genome string.
    fn get_invalid_characters(&'a self) -> Vec<u8> {
        self.iter()
            .copied()
            .filter(|c| !is_valid_ascii_genome_character(*c))
            .unique()
            .collect()
    }

    /// Copies this genome string into a `Vec`.
    fn clone_as_vec(&'a self) -> Vec<u8> {
        self.iter().copied().collect()
    }

    /// Get a reference to this genome as its subsequence type.
    fn as_genome_subsequence(&self) -> &GenomeSubsequence {
        println!("as_genome_subsequence()");
        self.index(0..self.len())
    }

    /// Returns the genome as nucleotide string.
    fn as_string(&'a self) -> String {
        String::from_utf8(self.clone_as_vec())
            .expect("Genome contains non-utf8 characters (It should be ASCII only).")
    }

    /// Returns an iterator over the reverse complement of this genome.
    /// Panics if the iterator his an invalid character (see [not valid](is_valid)).
    fn reverse_complement_iter(&'a self) -> ReverseComplementIterator<Self::Iterator> {
        self.iter()
            .copied()
            .rev()
            .map(ascii_complement as fn(u8) -> Option<u8>)
            .map(Option::unwrap as fn(Option<u8>) -> u8)
    }
}

/// A genome sequence that is owned, i.e. not a reference.
pub trait OwnedGenomeSequence<'a, GenomeSubsequence: GenomeSequence<'a, GenomeSubsequence> + ?Sized>:
    for<'s> GenomeSequence<'s, GenomeSubsequence> + FromIterator<u8>
{
    /// Returns the reverse complement of this genome.
    /// Panics if this genome is [not valid](is_valid).
    fn reverse_complement(&'a self) -> Self {
        self.reverse_complement_iter().collect()
    }
}

/// A mutable genome sequence.
pub trait GenomeSequenceMut<
    'a,
    GenomeSubsequenceMut: GenomeSequenceMut<'a, GenomeSubsequenceMut> + ?Sized,
>: SequenceMut<'a, u8, GenomeSubsequenceMut> + GenomeSequence<'a, GenomeSubsequenceMut>
{
}

/// An editable genome sequence.
pub trait EditableGenomeSequence<
    'a,
    GenomeSubsequence: GenomeSequence<'a, GenomeSubsequence> + ?Sized,
>: EditableSequence<'a, u8, GenomeSubsequence> + GenomeSequence<'a, GenomeSubsequence>
{
}

/// Returns the complement of the given genome char.
/// Returns `None` if the given char [is invalid](is_valid_ascii_genome_character).
pub fn ascii_complement(char: u8) -> Option<u8> {
    match char {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'G' => Some(b'C'),
        b'C' => Some(b'G'),
        _ => None,
    }
}

/// Returns true if the given ascii character represents a valid genome character.
/// Valid genome characters are `A`, `T`, `G` and `C`.
// Note: do not add more characters here, but make a new method if required.
pub fn is_valid_ascii_genome_character(char: u8) -> bool {
    matches!(char, b'A' | b'T' | b'G' | b'C')
}
