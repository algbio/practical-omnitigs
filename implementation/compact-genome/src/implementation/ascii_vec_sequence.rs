//! A simple representation of a genome as `Vec<u8>`.

use crate::interface::sequence::{
    EditableGenomeSequence, GenomeSequence, GenomeSequenceMut, OwnedGenomeSequence,
};

/// A genome sequence stored as vector of ASCII characters.
pub type AsciiVectorGenome = Vec<u8>;

impl<'a> GenomeSequence<'a, [u8]> for Vec<u8> {}

impl<'a> GenomeSequenceMut<'a, [u8]> for Vec<u8> {}

impl<'a> OwnedGenomeSequence<'a, [u8]> for Vec<u8> {}

impl<'a> EditableGenomeSequence<'a, [u8]> for Vec<u8> {}

impl<'a> GenomeSequence<'a, [u8]> for [u8] {}

impl<'a> GenomeSequenceMut<'a, [u8]> for [u8] {}

#[cfg(test)]
mod tests {
    use crate::implementation::ascii_vec_sequence::AsciiVectorGenome;
    use crate::interface::sequence::{GenomeSequence, OwnedGenomeSequence};

    #[test]
    fn test_reverse_complement() {
        let genome: AsciiVectorGenome = b"ATTCGGT".iter().copied().collect();
        let reverse_complement: AsciiVectorGenome = b"ACCGAAT".iter().copied().collect();
        debug_assert_eq!(genome.clone_as_reverse_complement(), reverse_complement);
        debug_assert_eq!(genome, reverse_complement.clone_as_reverse_complement());
    }

    #[test]
    fn test_display() {
        let genome: AsciiVectorGenome = b"ATTCGGT".iter().copied().collect();
        let display_string = genome.as_string();
        let expected_string = "ATTCGGT";
        debug_assert_eq!(display_string, expected_string);
    }
}
