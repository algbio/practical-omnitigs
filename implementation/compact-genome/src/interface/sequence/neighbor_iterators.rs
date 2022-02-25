//! Iterators over different neighborhoods of genome sequences.

use crate::interface::alphabet::{Alphabet, AlphabetCharacter};
use crate::interface::sequence::{GenomeSequenceMut, OwnedGenomeSequence};
use std::marker::PhantomData;

/// Returns an iterator over the genome sequences that are exactly one substitution away from this genome sequence.
/// In other words, all sequences with hamming distance one of the same length.
/// While the returned type works like an iterator, it does not implement the iterator trait due to limitations of the type system.
pub fn substitution_distance_one_neighbor_iterator<
    'owned_sequence,
    AlphabetType: Alphabet,
    Sequence: OwnedGenomeSequence<'owned_sequence, AlphabetType, Subsequence>
        + for<'sequence> GenomeSequenceMut<'sequence, AlphabetType, Subsequence>
        + Sized,
    Subsequence: 'owned_sequence
        + for<'subsequence> GenomeSequenceMut<'subsequence, AlphabetType, Subsequence>
        + ?Sized,
>(
    sequence: Sequence,
) -> SubstitutionDistanceOneNeighborIterator<'owned_sequence, AlphabetType, Sequence, Subsequence> {
    SubstitutionDistanceOneNeighborIterator::new(sequence)
}

/// An iterator over the genome sequences that are exactly one substitution away from the given genome sequence.
/// In other words, all sequences with hamming distance one of the same length.
/// While the type works like an iterator, it does not implement the iterator trait due to limitations of the type system.
pub struct SubstitutionDistanceOneNeighborIterator<
    'owned_sequence,
    AlphabetType: Alphabet,
    Sequence: OwnedGenomeSequence<'owned_sequence, AlphabetType, Subsequence>
        + for<'sequence> GenomeSequenceMut<'sequence, AlphabetType, Subsequence>
        + Sized,
    Subsequence: 'owned_sequence
        + for<'subsequence> GenomeSequenceMut<'subsequence, AlphabetType, Subsequence>
        + ?Sized,
> {
    current_index: usize,
    current_character: usize,
    original_character: AlphabetType::CharacterType,
    sequence: Sequence,
    subsequence: PhantomData<&'owned_sequence Subsequence>,
    alphabet_type: PhantomData<AlphabetType>,
}

impl<
        'owned_sequence,
        AlphabetType: Alphabet,
        Sequence: OwnedGenomeSequence<'owned_sequence, AlphabetType, Subsequence>
            + for<'sequence> GenomeSequenceMut<'sequence, AlphabetType, Subsequence>
            + Sized,
        Subsequence: 'owned_sequence
            + for<'subsequence> GenomeSequenceMut<'subsequence, AlphabetType, Subsequence>
            + ?Sized,
    >
    SubstitutionDistanceOneNeighborIterator<'owned_sequence, AlphabetType, Sequence, Subsequence>
{
    fn new(sequence: Sequence) -> Self {
        Self {
            current_index: 0,
            current_character: 0,
            original_character: if sequence.len() > 0 {
                sequence[0].clone()
            } else {
                AlphabetType::CharacterType::from_index(0).unwrap()
            },
            sequence,
            subsequence: Default::default(),
            alphabet_type: Default::default(),
        }
    }

    /// Like the `Iterator::next` function.
    /// The returned reference must be dropped before next is called a second time.
    pub fn next<'this: 'returned_reference, 'returned_reference>(
        &'this mut self,
    ) -> Option<&'returned_reference Subsequence> {
        while self.current_index < self.sequence.len() {
            while self.current_character < AlphabetType::SIZE {
                let current_character =
                    AlphabetType::CharacterType::from_index(self.current_character).unwrap();
                self.current_character += 1;

                if self.original_character != current_character {
                    self.sequence[self.current_index] = current_character;
                    return Some(self.sequence.as_genome_subsequence());
                }
            }

            self.sequence[self.current_index] = self.original_character.clone();
            self.current_index += 1;
            self.current_character = 0;
            if self.current_index < self.sequence.len() {
                self.original_character = self.sequence[self.current_index].clone();
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use crate::implementation::vec_sequence::VectorGenome;
    use crate::interface::alphabet::dna_alphabet::DnaAlphabet;
    use crate::interface::sequence::neighbor_iterators::substitution_distance_one_neighbor_iterator;
    use crate::interface::sequence::OwnedGenomeSequence;

    #[test]
    fn test_substitution_distance_one_neighbor_iterator() {
        let sequence = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ACCGTTA").unwrap();
        let mut neighbors = Vec::new();
        let mut neighbor_iterator = substitution_distance_one_neighbor_iterator(sequence);
        while let Some(neighbor) = neighbor_iterator.next() {
            neighbors.push(neighbor.to_owned());
        }
        neighbors.sort();
        let mut expected = vec![
            VectorGenome::from_slice_u8(b"CCCGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"GCCGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"TCCGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"AACGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"AGCGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"ATCGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACAGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACGGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACTGTTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCATTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCCTTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCTTTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGATA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGCTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGGTA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGTAA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGTCA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGTGA").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGTTC").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGTTG").unwrap(),
            VectorGenome::from_slice_u8(b"ACCGTTT").unwrap(),
        ];
        expected.sort();
        debug_assert_eq!(neighbors, expected);
    }

    #[test]
    fn test_substitution_distance_one_neighbor_iterator_empty() {
        let sequence = VectorGenome::from_slice_u8(b"").unwrap();
        let mut neighbors = Vec::new();
        let mut neighbor_iterator = substitution_distance_one_neighbor_iterator(sequence);
        while let Some(neighbor) = neighbor_iterator.next() {
            neighbors.push(neighbor.to_owned());
        }
        neighbors.sort();
        let mut expected: Vec<VectorGenome<DnaAlphabet>> = Vec::new();
        expected.sort();
        debug_assert_eq!(neighbors, expected);
    }
}
