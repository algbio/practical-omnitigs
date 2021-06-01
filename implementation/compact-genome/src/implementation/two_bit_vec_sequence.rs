//! A representation of a genome as `Vec<usize>` where each character is encoded as two bits.

use crate::interface::sequence::{EditableGenomeSequence, GenomeSequence, OwnedGenomeSequence};
use bitvec::prelude::*;
use ref_cast::RefCast;
use std::borrow::Borrow;
use std::iter::FromIterator;
use std::ops::{Index, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive};
use traitsequence::interface::{EditableSequence, OwnedSequence, Sequence};

const ASCII_A: u8 = b'A';
const ASCII_C: u8 = b'C';
const ASCII_G: u8 = b'G';
const ASCII_T: u8 = b'T';

/// A genome sequence stored as vector of two-bit characters.
///
/// Translation table:
/// ```txt
/// 00 - A
/// 01 - C
/// 10 - G
/// 11 - T
/// ```
#[derive(Default, Debug, Clone, Eq, PartialEq, Hash)]
pub struct TwoBitVectorGenome {
    /// Stores the sequence as two-bit characters.
    bits: BitVec,
}

/// The subsequence of a genome sequence stored as vector of two-bit characters.
#[derive(RefCast, Debug, Eq, PartialEq, Hash)]
#[repr(transparent)]
pub struct TwoBitVectorSubGenome {
    bits: BitSlice,
}

/// An iterator over a [TwoBitVectorGenome].
pub struct TwoBitVectorGenomeIterator {
    current: usize,
    sequence: TwoBitVectorGenome,
}

/// An iterator over a [TwoBitVectorSubGenome].
pub struct TwoBitVectorSubGenomeIterator<'a> {
    slice: &'a TwoBitVectorSubGenome,
}

impl TwoBitVectorGenome {
    /// Returns the amount of memory this genome sequence uses in bytes.
    /// This is meant to be accurate, but might be off by a constant number of bytes.
    pub fn size_in_memory(&self) -> usize {
        std::mem::size_of::<BitVec>() + self.bits.capacity() / 8
    }
}

impl<'a> Sequence<'a, u8, TwoBitVectorSubGenome> for TwoBitVectorGenome {
    type Iterator = TwoBitVectorSubGenomeIterator<'a>;

    fn iter(&'a self) -> Self::Iterator {
        self.as_genome_subsequence().iter()
    }

    fn len(&self) -> usize {
        self.as_genome_subsequence().len()
    }
}

impl<'a> Sequence<'a, u8, TwoBitVectorSubGenome> for TwoBitVectorSubGenome {
    type Iterator = TwoBitVectorSubGenomeIterator<'a>;

    fn iter(&'a self) -> Self::Iterator {
        TwoBitVectorSubGenomeIterator { slice: self }
    }

    fn len(&self) -> usize {
        self.bits.len() / 2
    }
}

impl<'a> EditableSequence<'a, u8, TwoBitVectorSubGenome> for TwoBitVectorGenome {}

impl Index<Range<usize>> for TwoBitVectorGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl Index<RangeFrom<usize>> for TwoBitVectorGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeFrom<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl Index<RangeTo<usize>> for TwoBitVectorGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeTo<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl Index<RangeFull> for TwoBitVectorGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeFull) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl Index<RangeInclusive<usize>> for TwoBitVectorGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeInclusive<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl Index<RangeToInclusive<usize>> for TwoBitVectorGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeToInclusive<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl Index<usize> for TwoBitVectorGenome {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl Index<Range<usize>> for TwoBitVectorSubGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        TwoBitVectorSubGenome::ref_cast(&self.bits[index.start * 2..index.end * 2])
    }
}

impl Index<RangeFrom<usize>> for TwoBitVectorSubGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeFrom<usize>) -> &Self::Output {
        self.index(index.start..self.len())
    }
}

impl Index<RangeTo<usize>> for TwoBitVectorSubGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeTo<usize>) -> &Self::Output {
        self.index(0..index.end)
    }
}

impl Index<RangeFull> for TwoBitVectorSubGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, _index: RangeFull) -> &Self::Output {
        self.index(0..self.len())
    }
}

impl Index<RangeInclusive<usize>> for TwoBitVectorSubGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeInclusive<usize>) -> &Self::Output {
        TwoBitVectorSubGenome::ref_cast(&self.bits[index.start() * 2..=index.end() * 2])
    }
}

impl Index<RangeToInclusive<usize>> for TwoBitVectorSubGenome {
    type Output = TwoBitVectorSubGenome;

    fn index(&self, index: RangeToInclusive<usize>) -> &Self::Output {
        self.index(0..=index.end)
    }
}

impl Index<usize> for TwoBitVectorSubGenome {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        match (self.bits[index * 2 + 1], self.bits[index * 2]) {
            (false, false) => &ASCII_A,
            (false, true) => &ASCII_C,
            (true, false) => &ASCII_G,
            (true, true) => &ASCII_T,
        }
    }
}

impl FromIterator<u8> for TwoBitVectorGenome {
    fn from_iter<T: IntoIterator<Item = u8>>(iter: T) -> Self {
        let mut result = Self::default();
        result.extend(iter);
        result
    }
}

impl Extend<u8> for TwoBitVectorGenome {
    fn extend<T: IntoIterator<Item = u8>>(&mut self, iter: T) {
        for ascii_char in iter {
            match ascii_char {
                ASCII_A => {
                    self.bits.push(false);
                    self.bits.push(false);
                }
                ASCII_C => {
                    self.bits.push(true);
                    self.bits.push(false);
                }
                ASCII_G => {
                    self.bits.push(false);
                    self.bits.push(true);
                }
                ASCII_T => {
                    self.bits.push(true);
                    self.bits.push(true);
                }
                c => panic!("Invalid character {} with ASCII code {}", c as char, c),
            }
        }
    }
}

impl IntoIterator for TwoBitVectorGenome {
    type Item = u8;
    type IntoIter = TwoBitVectorGenomeIterator;

    fn into_iter(self) -> Self::IntoIter {
        TwoBitVectorGenomeIterator {
            sequence: self,
            current: 0,
        }
    }
}

impl Iterator for TwoBitVectorGenomeIterator {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.sequence.len() {
            let result = self.sequence[self.current];
            self.current += 1;
            Some(result)
        } else {
            None
        }
    }
}

impl<'a> Iterator for TwoBitVectorSubGenomeIterator<'a> {
    type Item = &'a u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.slice.len() > 0 {
            let result = &self.slice[0];
            self.slice = &self.slice[1..self.slice.len()];
            Some(result)
        } else {
            None
        }
    }
}

impl<'a> DoubleEndedIterator for TwoBitVectorSubGenomeIterator<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.slice.len() > 0 {
            let result = &self.slice[self.slice.len() - 1];
            self.slice = &self.slice[0..self.slice.len() - 1];
            Some(result)
        } else {
            None
        }
    }
}

impl Borrow<TwoBitVectorSubGenome> for TwoBitVectorGenome {
    fn borrow(&self) -> &TwoBitVectorSubGenome {
        self.as_genome_subsequence()
    }
}

impl ToOwned for TwoBitVectorSubGenome {
    type Owned = TwoBitVectorGenome;

    fn to_owned(&self) -> Self::Owned {
        self.iter().copied().collect()
    }
}

impl<'a> GenomeSequence<'a, TwoBitVectorSubGenome> for TwoBitVectorGenome {
    fn as_genome_subsequence(&self) -> &TwoBitVectorSubGenome {
        TwoBitVectorSubGenome::ref_cast(&self.bits[..])
    }
}

impl<'a> OwnedGenomeSequence<'a, TwoBitVectorSubGenome> for TwoBitVectorGenome {}

impl<'a> EditableGenomeSequence<'a, TwoBitVectorSubGenome> for TwoBitVectorGenome {}

impl<'a> GenomeSequence<'a, TwoBitVectorSubGenome> for TwoBitVectorSubGenome {}

impl<'a> OwnedSequence<'a, u8, TwoBitVectorSubGenome> for TwoBitVectorGenome {}

#[cfg(test)]
mod tests {
    use crate::implementation::two_bit_vec_sequence::TwoBitVectorGenome;
    use crate::interface::sequence::{GenomeSequence, OwnedGenomeSequence};

    #[test]
    fn test_reverse_complement() {
        let genome: TwoBitVectorGenome = b"ATTCGGT".iter().copied().collect();
        let reverse_complement: TwoBitVectorGenome = b"ACCGAAT".iter().copied().collect();
        assert_eq!(genome.clone_as_reverse_complement(), reverse_complement);
        assert_eq!(genome, reverse_complement.clone_as_reverse_complement());
    }

    #[test]
    fn test_display() {
        let genome: TwoBitVectorGenome = b"ATTCGGT".iter().copied().collect();
        let display_string = genome.as_string();
        let expected_string = "ATTCGGT";
        assert_eq!(display_string, expected_string);
    }

    #[test]
    fn test_substrings() {
        let genome: TwoBitVectorGenome = b"ATTCGGT".iter().copied().collect();

        let display_string = genome[1..4].as_string();
        let expected_string = "TTC";
        assert_eq!(display_string, expected_string);

        let display_string = genome.clone_as_reverse_complement()[1..4].as_string();
        let expected_string = "CCG";
        assert_eq!(display_string, expected_string);

        let display_string =
            genome[1..6].to_owned().clone_as_reverse_complement()[1..4].as_string();
        let expected_string = "CGA";
        assert_eq!(display_string, expected_string);
    }

    #[test]
    fn test_empty_substring_after_end() {
        let genome: TwoBitVectorGenome = b"ATTCGGT".iter().copied().collect();
        let display_string = genome[7..7].as_string();
        let expected_string = "";
        assert_eq!(display_string, expected_string);
    }

    #[test]
    #[should_panic]
    fn test_empty_substring_after_end2() {
        let genome: TwoBitVectorGenome = b"ATTCGGT".iter().copied().collect();
        let display_string = genome[8..8].as_string();
        let expected_string = "";
        assert_eq!(display_string, expected_string);
    }
}
