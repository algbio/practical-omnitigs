//! A simple representation of a genome as `Vec<u8>`.

use crate::interface::{GenomeSequenceMut, GenomeSequence, OwnedGenomeSequence, EditableGenomeSequence};

/// A genome sequence stored as vector of ASCII characters.
pub type AsciiVectorGenome = Vec<u8>;

impl<'a> GenomeSequence<'a, [u8]> for Vec<u8> {
}

impl<'a> GenomeSequenceMut<'a, [u8]> for Vec<u8> {
}

impl<'a> OwnedGenomeSequence<'a, [u8]> for Vec<u8> {
}

impl<'a> EditableGenomeSequence<'a, [u8]> for Vec<u8> {
}

#[cfg(test)]
mod tests {
    use std::iter::FromIterator;
    use crate::interface::{OwnedGenomeSequence, GenomeSequence};

    #[test]
    fn test_reverse_complement() {
        let genome = Vec::from_iter(b"ATTCGGT".iter().copied());
        let reverse_complement = Vec::from_iter(b"ACCGAAT".iter().copied());
        assert_eq!(genome.reverse_complement(), reverse_complement);
        assert_eq!(genome, reverse_complement.reverse_complement());
    }

    #[test]
    fn test_display() {
        let genome = Vec::from_iter(b"ATTCGGT".iter().copied());
        let display_string = genome.as_string();
        let expected_string = "ATTCGGT";
        assert_eq!(display_string, expected_string);
    }
}