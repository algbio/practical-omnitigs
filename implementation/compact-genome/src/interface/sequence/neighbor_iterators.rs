//! Iterators over different neighborhoods of genome sequences.

use crate::interface::sequence::{GenomeSequenceMut, OwnedGenomeSequence};
use crate::{ASCII_A, ASCII_C, ASCII_G, ASCII_T};
use std::marker::PhantomData;

/// Returns an iterator over the genome sequences that are exactly one substitution away from this genome sequence.
/// In other words, all sequences with hamming distance one of the same length.
/// While the returned type works like an iterator, it does not implement the iterator trait due to limitations of the type system.
pub fn substitution_distance_one_neighbor_iterator<
    'owned_sequence,
    Sequence: OwnedGenomeSequence<'owned_sequence, Subsequence>
        + for<'sequence> GenomeSequenceMut<'sequence, Subsequence>
        + Sized,
    Subsequence: 'owned_sequence + for<'subsequence> GenomeSequenceMut<'subsequence, Subsequence> + ?Sized,
>(
    sequence: Sequence,
) -> SubstitutionDistanceOneNeighborIterator<'owned_sequence, Sequence, Subsequence> {
    SubstitutionDistanceOneNeighborIterator::new(sequence)
}

/// An iterator over the genome sequences that are exactly one substitution away from the given genome sequence.
/// In other words, all sequences with hamming distance one of the same length.
/// While the type works like an iterator, it does not implement the iterator trait due to limitations of the type system.
pub struct SubstitutionDistanceOneNeighborIterator<
    'owned_sequence,
    Sequence: OwnedGenomeSequence<'owned_sequence, Subsequence>
        + for<'sequence> GenomeSequenceMut<'sequence, Subsequence>
        + Sized,
    Subsequence: 'owned_sequence + for<'subsequence> GenomeSequenceMut<'subsequence, Subsequence> + ?Sized,
> {
    current_index: usize,
    current_character: u8,
    original_character: u8,
    sequence: Sequence,
    subsequence: PhantomData<&'owned_sequence Subsequence>,
}

impl<
        'owned_sequence,
        Sequence: OwnedGenomeSequence<'owned_sequence, Subsequence>
            + for<'sequence> GenomeSequenceMut<'sequence, Subsequence>
            + Sized,
        Subsequence: 'owned_sequence + for<'subsequence> GenomeSequenceMut<'subsequence, Subsequence> + ?Sized,
    > SubstitutionDistanceOneNeighborIterator<'owned_sequence, Sequence, Subsequence>
{
    fn new(sequence: Sequence) -> Self {
        Self {
            current_index: 0,
            current_character: 0,
            original_character: if sequence.len() > 0 { sequence[0] } else { 0 },
            sequence,
            subsequence: Default::default(),
        }
    }

    /// Like the `Iterator::next` function.
    /// The returned reference must be dropped before next is called a second time.
    pub fn next<'this: 'returned_reference, 'returned_reference>(
        &'this mut self,
    ) -> Option<&'returned_reference Subsequence> {
        while self.current_index < self.sequence.len() {
            while self.current_character < 4 {
                let current_ascii = match self.current_character {
                    0 => ASCII_A,
                    1 => ASCII_C,
                    2 => ASCII_G,
                    3 => ASCII_T,
                    _ => unreachable!(),
                };
                self.current_character += 1;

                if self.original_character != current_ascii {
                    self.sequence[self.current_index] = current_ascii;
                    return Some(self.sequence.as_genome_subsequence());
                }
            }

            self.sequence[self.current_index] = self.original_character;
            self.current_index += 1;
            self.current_character = 0;
            if self.current_index < self.sequence.len() {
                self.original_character = self.sequence[self.current_index];
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use crate::interface::sequence::neighbor_iterators::substitution_distance_one_neighbor_iterator;

    #[test]
    fn test_substitution_distance_one_neighbor_iterator() {
        let sequence = b"ACCGTTA".to_vec();
        let mut neighbors = Vec::new();
        let mut neighbor_iterator = substitution_distance_one_neighbor_iterator(sequence);
        while let Some(neighbor) = neighbor_iterator.next() {
            neighbors.push(neighbor.to_owned());
        }
        neighbors.sort();
        let mut expected = vec![
            b"CCCGTTA".to_vec(),
            b"GCCGTTA".to_vec(),
            b"TCCGTTA".to_vec(),
            b"AACGTTA".to_vec(),
            b"AGCGTTA".to_vec(),
            b"ATCGTTA".to_vec(),
            b"ACAGTTA".to_vec(),
            b"ACGGTTA".to_vec(),
            b"ACTGTTA".to_vec(),
            b"ACCATTA".to_vec(),
            b"ACCCTTA".to_vec(),
            b"ACCTTTA".to_vec(),
            b"ACCGATA".to_vec(),
            b"ACCGCTA".to_vec(),
            b"ACCGGTA".to_vec(),
            b"ACCGTAA".to_vec(),
            b"ACCGTCA".to_vec(),
            b"ACCGTGA".to_vec(),
            b"ACCGTTC".to_vec(),
            b"ACCGTTG".to_vec(),
            b"ACCGTTT".to_vec(),
        ];
        expected.sort();
        debug_assert_eq!(neighbors, expected);
    }

    #[test]
    fn test_substitution_distance_one_neighbor_iterator_empty() {
        let sequence = b"".to_vec();
        let mut neighbors = Vec::new();
        let mut neighbor_iterator = substitution_distance_one_neighbor_iterator(sequence);
        while let Some(neighbor) = neighbor_iterator.next() {
            neighbors.push(neighbor.to_owned());
        }
        neighbors.sort();
        let mut expected: Vec<Vec<u8>> = Vec::new();
        expected.sort();
        debug_assert_eq!(neighbors, expected);
    }
}
