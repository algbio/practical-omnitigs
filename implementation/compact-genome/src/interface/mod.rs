//! The traits providing the abstractions of this crate.

use itertools::Itertools;
use std::iter::FromIterator;
use traitsequence::interface::Sequence;

/// A genome string.
/// While the internal representation is implementation specific, externally genome strings are represented as sequences of `u8`,
/// which ASCII encode the valid genome characters specified by [is_valid_ascii_genome_character()](is_valid_ascii_genome_character).
///
/// The ordering implemented should be the lexical order of the genome characters.
pub trait Genome:
    for<'a> FromIterator<&'a u8>
    + FromIterator<u8>
    + Eq
    + Clone
    + Ord
    + Sized
    + for<'a> Sequence<'a, u8>
where
    for<'a> &'a Self: IntoIterator<Item = u8>,
{
    /// Returns the reverse complement of this genome.
    /// Panics if this genome is [not valid](is_valid).
    fn reverse_complement(&self) -> Self;

    /// Returns true if this genome is valid, i.e. it contains no invalid characters.
    /// Valid characters are defined by [is_valid_ascii_genome_character()](is_valid_ascii_genome_character)
    fn is_valid(&self) -> bool {
        self.into_iter().all(is_valid_ascii_genome_character)
    }

    /// Returns a duplicate-free vector of all invalid characters in this genome string.
    fn get_invalid_characters(&self) -> Vec<u8> {
        self.into_iter()
            .filter(|c| !is_valid_ascii_genome_character(*c))
            .unique()
            .collect()
    }

    /// Copies this genome string into a `Vec`.
    fn clone_as_vec(&self) -> Vec<u8> {
        self.into_iter().collect()
    }

    /// Returns a copy of the prefix with length `len` of this genome.
    /// Panics if `len >= self.len()`.
    fn prefix(&self, len: usize) -> Self {
        Self::from_iter(self.into_iter().take(len))
    }

    /// Returns a copy of the suffix with length `len` of this genome.
    /// Panics if `len >= self.len()`.
    fn suffix(&self, len: usize) -> Self {
        Self::from_iter(self.into_iter().skip(self.len() - len))
    }

    /// Returns the genome as nucleotide string.
    fn as_string(&self) -> String {
        String::from_utf8(self.clone_as_vec())
            .expect("Genome contains non-utf8 characters (It should be ASCII only).")
    }
}

/// A genome string that can be extended with another genome string.
pub trait ExtendableGenome: Genome
where
    for<'a> &'a Self: IntoIterator<Item = u8>,
{
    /// Append the elements of the given source to this genome.
    fn extend<ExtensionSource: IntoIterator<Item = u8>>(&mut self, extension: ExtensionSource);
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
