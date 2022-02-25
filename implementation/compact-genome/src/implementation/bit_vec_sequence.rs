//! A representation of a genome as `Vec<usize>` where each character is encoded as bits.

use crate::interface::alphabet::{Alphabet, AlphabetCharacter};
use crate::interface::sequence::{EditableGenomeSequence, GenomeSequence, OwnedGenomeSequence};
use bitvec::prelude::*;
use ref_cast::RefCast;
use std::borrow::Borrow;
use std::iter::FromIterator;
use std::marker::PhantomData;
use std::mem;
use std::ops::{Index, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive};
use traitsequence::interface::{EditableSequence, OwnedSequence, Sequence};

/// A genome sequence stored as vector of minimum-bit characters.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct BitVectorGenome<AlphabetType: Alphabet> {
    phantom_data: PhantomData<AlphabetType>,
    /// Stores the sequence as minimum-bit characters.
    pub(crate) bits: BitVec,
}

/// The subsequence of a genome sequence stored as vector of minimum-bit characters.
#[derive(RefCast, Debug, Eq, PartialEq, Hash)]
#[repr(transparent)]
pub struct BitVectorSubGenome<AlphabetType: Alphabet> {
    phantom_data: PhantomData<AlphabetType>,
    pub(crate) bits: BitSlice,
}

/// An iterator over a [BitVectorGenome].
pub struct BitVectorGenomeIterator<AlphabetType: Alphabet> {
    current: usize,
    sequence: BitVectorGenome<AlphabetType>,
}

/// An iterator over a [BitVectorSubGenome].
pub struct BitVectorSubGenomeIterator<'a, AlphabetType: Alphabet> {
    slice: &'a BitVectorSubGenome<AlphabetType>,
}

impl<AlphabetType: Alphabet> BitVectorGenome<AlphabetType> {
    /// Returns the amount of memory this genome sequence uses in bytes.
    /// This is meant to be accurate, but might be off by a constant number of bytes.
    pub fn size_in_memory(&self) -> usize {
        std::mem::size_of::<BitVec>() + self.bits.capacity() / 8
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    Sequence<'a, AlphabetType::CharacterType, BitVectorSubGenome<AlphabetType>>
    for BitVectorGenome<AlphabetType>
{
    type Iterator = BitVectorSubGenomeIterator<'a, AlphabetType>;

    fn iter(&'a self) -> Self::Iterator {
        self.as_genome_subsequence().iter()
    }

    fn len(&self) -> usize {
        self.as_genome_subsequence().len()
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    Sequence<'a, AlphabetType::CharacterType, BitVectorSubGenome<AlphabetType>>
    for BitVectorSubGenome<AlphabetType>
{
    type Iterator = BitVectorSubGenomeIterator<'a, AlphabetType>;

    fn iter(&'a self) -> Self::Iterator {
        BitVectorSubGenomeIterator { slice: self }
    }

    fn len(&self) -> usize {
        self.bits.len() / alphabet_character_bit_width(AlphabetType::SIZE)
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    EditableSequence<'a, AlphabetType::CharacterType, BitVectorSubGenome<AlphabetType>>
    for BitVectorGenome<AlphabetType>
{
}

impl<AlphabetType: Alphabet> Index<Range<usize>> for BitVectorGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeFrom<usize>> for BitVectorGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeFrom<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeTo<usize>> for BitVectorGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeTo<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeFull> for BitVectorGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeFull) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeInclusive<usize>> for BitVectorGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeInclusive<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeToInclusive<usize>> for BitVectorGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeToInclusive<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<usize> for BitVectorGenome<AlphabetType> {
    type Output = AlphabetType::CharacterType;

    fn index(&self, index: usize) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<Range<usize>> for BitVectorSubGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        BitVectorSubGenome::ref_cast(&self.bits[index.start * bit_width..index.end * bit_width])
    }
}

impl<AlphabetType: Alphabet> Index<RangeFrom<usize>> for BitVectorSubGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeFrom<usize>) -> &Self::Output {
        self.index(index.start..self.len())
    }
}

impl<AlphabetType: Alphabet> Index<RangeTo<usize>> for BitVectorSubGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeTo<usize>) -> &Self::Output {
        self.index(0..index.end)
    }
}

impl<AlphabetType: Alphabet> Index<RangeFull> for BitVectorSubGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, _index: RangeFull) -> &Self::Output {
        self.index(0..self.len())
    }
}

impl<AlphabetType: Alphabet> Index<RangeInclusive<usize>> for BitVectorSubGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeInclusive<usize>) -> &Self::Output {
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        BitVectorSubGenome::ref_cast(
            &self.bits[index.start() * bit_width..=index.end() * bit_width],
        )
    }
}

impl<AlphabetType: Alphabet> Index<RangeToInclusive<usize>> for BitVectorSubGenome<AlphabetType> {
    type Output = BitVectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeToInclusive<usize>) -> &Self::Output {
        self.index(0..=index.end)
    }
}

impl<AlphabetType: Alphabet> Index<usize> for BitVectorSubGenome<AlphabetType> {
    type Output = AlphabetType::CharacterType;

    fn index(&self, index: usize) -> &Self::Output {
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        let offset = index * bit_width;
        let limit = (index + 1) * bit_width;
        let value: usize = self.bits[offset..limit].load();
        Self::Output::from_index_ref(value).expect("bitvec contains invalid character")
    }
}

impl<AlphabetType: Alphabet> FromIterator<AlphabetType::CharacterType>
    for BitVectorGenome<AlphabetType>
{
    fn from_iter<T: IntoIterator<Item = AlphabetType::CharacterType>>(iter: T) -> Self {
        let mut result = Self::default();
        result.extend(iter);
        result
    }
}

impl<AlphabetType: Alphabet> Extend<AlphabetType::CharacterType> for BitVectorGenome<AlphabetType> {
    fn extend<T: IntoIterator<Item = AlphabetType::CharacterType>>(&mut self, iter: T) {
        let iter = iter.into_iter();
        let (size, _) = iter.size_hint();
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        self.bits.reserve(size * bit_width);

        for character in iter {
            let value = character.index();
            self.bits
                .extend_from_bitslice(&value.view_bits::<Lsb0>()[0..bit_width]);
        }
    }
}

impl<AlphabetType: Alphabet> IntoIterator for BitVectorGenome<AlphabetType> {
    type Item = AlphabetType::CharacterType;
    type IntoIter = BitVectorGenomeIterator<AlphabetType>;

    fn into_iter(self) -> Self::IntoIter {
        BitVectorGenomeIterator {
            sequence: self,
            current: 0,
        }
    }
}

impl<AlphabetType: Alphabet> Iterator for BitVectorGenomeIterator<AlphabetType> {
    type Item = AlphabetType::CharacterType;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.sequence.len() {
            let result = self.sequence[self.current].clone();
            self.current += 1;
            Some(result)
        } else {
            None
        }
    }
}

impl<'iter, AlphabetType: Alphabet> Iterator for BitVectorSubGenomeIterator<'iter, AlphabetType> {
    type Item = &'iter AlphabetType::CharacterType;

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

impl<'a, AlphabetType: Alphabet> DoubleEndedIterator
    for BitVectorSubGenomeIterator<'a, AlphabetType>
{
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

impl<AlphabetType: Alphabet> Borrow<BitVectorSubGenome<AlphabetType>>
    for BitVectorGenome<AlphabetType>
{
    fn borrow(&self) -> &BitVectorSubGenome<AlphabetType> {
        self.as_genome_subsequence()
    }
}

impl<AlphabetType: Alphabet> ToOwned for BitVectorSubGenome<AlphabetType> {
    type Owned = BitVectorGenome<AlphabetType>;

    fn to_owned(&self) -> Self::Owned {
        self.iter().cloned().collect()
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    GenomeSequence<'a, AlphabetType, BitVectorSubGenome<AlphabetType>>
    for BitVectorGenome<AlphabetType>
{
    fn as_genome_subsequence(&self) -> &BitVectorSubGenome<AlphabetType> {
        BitVectorSubGenome::ref_cast(&self.bits[..])
    }
}

impl<'a, AlphabetType: Alphabet + 'static>
    OwnedGenomeSequence<'a, AlphabetType, BitVectorSubGenome<AlphabetType>>
    for BitVectorGenome<AlphabetType>
{
}

impl<'a, AlphabetType: Alphabet + 'a>
    EditableGenomeSequence<'a, AlphabetType, BitVectorSubGenome<AlphabetType>>
    for BitVectorGenome<AlphabetType>
{
    fn reserve(&mut self, additional: usize) {
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        self.bits.reserve(additional * bit_width)
    }

    fn resize(&mut self, len: usize, default: AlphabetType::CharacterType) {
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        if self.len() <= len {
            self.bits.resize(len * bit_width, false);
        } else {
            let difference = len - self.len();
            let value = default.index();
            for _ in 0..difference {
                self.bits
                    .extend_from_bitslice(&value.view_bits::<Lsb0>()[0..bit_width]);
            }
        }
    }

    fn push(&mut self, character: AlphabetType::CharacterType) {
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        let value = character.index();
        self.bits
            .extend_from_bitslice(&value.view_bits::<Lsb0>()[0..bit_width])
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    GenomeSequence<'a, AlphabetType, BitVectorSubGenome<AlphabetType>>
    for BitVectorSubGenome<AlphabetType>
{
}

impl<'a, AlphabetType: Alphabet + 'a>
    OwnedSequence<'a, AlphabetType::CharacterType, BitVectorSubGenome<AlphabetType>>
    for BitVectorGenome<AlphabetType>
{
}

impl<AlphabetType: Alphabet> Default for BitVectorGenome<AlphabetType> {
    fn default() -> Self {
        Self {
            phantom_data: Default::default(),
            bits: Default::default(),
        }
    }
}

pub(crate) const fn alphabet_character_bit_width(size: usize) -> usize {
    mem::size_of::<usize>() * 8 - ((size - 1).leading_zeros() as usize)
}

#[cfg(test)]
mod tests {
    use crate::implementation::bit_vec_sequence::BitVectorGenome;
    use crate::interface::alphabet::dna_alphabet::DnaAlphabet;
    use crate::interface::alphabet::Alphabet;
    use crate::interface::sequence::{EditableGenomeSequence, GenomeSequence, OwnedGenomeSequence};

    #[test]
    fn test_reverse_complement() {
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();
        let reverse_complement = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ACCGAAT").unwrap();
        assert_eq!(genome.clone_as_reverse_complement(), reverse_complement);
        assert_eq!(genome, reverse_complement.clone_as_reverse_complement());
    }

    #[test]
    fn test_display() {
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();
        let display_string = genome.as_string();
        let expected_string = "ATTCGGT";
        assert_eq!(display_string, expected_string);
    }

    #[test]
    fn test_substrings() {
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();

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
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();
        let display_string = genome[7..7].as_string();
        let expected_string = "";
        assert_eq!(display_string, expected_string);
    }

    #[test]
    #[should_panic]
    fn test_empty_substring_after_end2() {
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();
        let display_string = genome[8..8].as_string();
        let expected_string = "";
        assert_eq!(display_string, expected_string);
    }

    #[test]
    fn test_canonical() {
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();
        assert!(!genome.is_canonical());
        assert!(genome.clone_as_reverse_complement().is_canonical());
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATAT").unwrap();
        assert!(genome.is_canonical());
        assert!(genome.clone_as_reverse_complement().is_canonical());
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"CGTA").unwrap();
        assert!(genome.is_canonical());
        assert!(!genome.clone_as_reverse_complement().is_canonical());
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"").unwrap();
        assert!(genome.is_canonical());
        assert!(genome.clone_as_reverse_complement().is_canonical());
    }

    #[test]
    fn test_self_complemental() {
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();
        assert!(!genome.is_self_complemental());
        assert!(!genome.clone_as_reverse_complement().is_self_complemental());
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATAT").unwrap();
        assert!(genome.is_self_complemental());
        assert!(genome.clone_as_reverse_complement().is_self_complemental());
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"CGTA").unwrap();
        assert!(!genome.is_self_complemental());
        assert!(!genome.clone_as_reverse_complement().is_self_complemental());
        let genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"").unwrap();
        assert!(genome.is_self_complemental());
        assert!(genome.clone_as_reverse_complement().is_self_complemental());
    }

    #[test]
    fn test_extend() {
        let mut genome = BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATT").unwrap();
        genome.extend_from_slice_u8(b"CGT").unwrap();
        assert_eq!(
            genome,
            BitVectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGT").unwrap()
        );
        assert_eq!(genome[0], DnaAlphabet::ascii_to_character(b'A').unwrap());
        assert_eq!(genome[1], DnaAlphabet::ascii_to_character(b'T').unwrap());
        assert_eq!(genome[2], DnaAlphabet::ascii_to_character(b'T').unwrap());
        assert_eq!(genome[3], DnaAlphabet::ascii_to_character(b'C').unwrap());
        assert_eq!(genome[4], DnaAlphabet::ascii_to_character(b'G').unwrap());
        assert_eq!(genome[5], DnaAlphabet::ascii_to_character(b'T').unwrap());
    }
}
