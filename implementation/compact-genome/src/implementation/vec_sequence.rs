//! A simple representation of a genome as `Vec<u8>`.

use std::borrow::Borrow;
use std::ops::{
    Index, IndexMut, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive,
};
//use std::marker::PhantomData;
use crate::interface::alphabet::Alphabet;
use crate::interface::sequence::{
    EditableGenomeSequence, GenomeSequence, GenomeSequenceMut, OwnedGenomeSequence,
};
use ref_cast::RefCast;
use traitsequence::interface::{EditableSequence, OwnedSequence, Sequence, SequenceMut};

/// A genome sequence stored as vector of plain characters.
#[derive(Debug, Clone, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct VectorGenome<AlphabetType: Alphabet> {
    vector: Vec<AlphabetType::CharacterType>,
}

/// The subsequence of a genome sequence stored as vector of plain characters.
#[derive(RefCast, Debug, Eq, PartialEq, Hash)]
#[repr(transparent)]
pub struct VectorSubGenome<AlphabetType: Alphabet> {
    //phantom_data: PhantomData<AlphabetType>,
    pub(crate) slice: [AlphabetType::CharacterType],
}

impl<'a, AlphabetType: Alphabet + 'a>
    GenomeSequence<'a, AlphabetType, VectorSubGenome<AlphabetType>> for VectorGenome<AlphabetType>
{
    fn as_genome_subsequence(&self) -> &VectorSubGenome<AlphabetType> {
        VectorSubGenome::ref_cast(&self.vector[..])
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    GenomeSequenceMut<'a, AlphabetType, VectorSubGenome<AlphabetType>>
    for VectorGenome<AlphabetType>
{
}

impl<'a, AlphabetType: Alphabet + 'static>
    OwnedGenomeSequence<'a, AlphabetType, VectorSubGenome<AlphabetType>>
    for VectorGenome<AlphabetType>
{
}

impl<'a, AlphabetType: Alphabet + 'a>
    EditableGenomeSequence<'a, AlphabetType, VectorSubGenome<AlphabetType>>
    for VectorGenome<AlphabetType>
{
    fn reserve(&mut self, additional: usize) {
        self.vector.reserve(additional)
    }

    fn resize(&mut self, len: usize, default: AlphabetType::CharacterType) {
        self.vector.resize(len, default)
    }

    fn push(&mut self, character: AlphabetType::CharacterType) {
        self.vector.push(character)
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    GenomeSequence<'a, AlphabetType, VectorSubGenome<AlphabetType>>
    for VectorSubGenome<AlphabetType>
{
}

impl<'a, AlphabetType: Alphabet + 'a>
    GenomeSequenceMut<'a, AlphabetType, VectorSubGenome<AlphabetType>>
    for VectorSubGenome<AlphabetType>
{
}

impl<'a, AlphabetType: Alphabet + 'a>
    Sequence<'a, AlphabetType::CharacterType, VectorSubGenome<AlphabetType>>
    for VectorGenome<AlphabetType>
{
    type Iterator = std::slice::Iter<'a, AlphabetType::CharacterType>;

    fn iter(&'a self) -> Self::Iterator {
        self.as_genome_subsequence().iter()
    }

    fn len(&self) -> usize {
        self.as_genome_subsequence().len()
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    Sequence<'a, AlphabetType::CharacterType, VectorSubGenome<AlphabetType>>
    for VectorSubGenome<AlphabetType>
{
    type Iterator = std::slice::Iter<'a, AlphabetType::CharacterType>;

    fn iter(&'a self) -> Self::Iterator {
        self.slice.iter()
    }

    fn len(&self) -> usize {
        self.slice.len()
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    EditableSequence<'a, AlphabetType::CharacterType, VectorSubGenome<AlphabetType>>
    for VectorGenome<AlphabetType>
{
}

impl<AlphabetType: Alphabet> Index<Range<usize>> for VectorGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeFrom<usize>> for VectorGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeFrom<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeTo<usize>> for VectorGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeTo<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeFull> for VectorGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeFull) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeInclusive<usize>> for VectorGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeInclusive<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<RangeToInclusive<usize>> for VectorGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeToInclusive<usize>) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<usize> for VectorGenome<AlphabetType> {
    type Output = AlphabetType::CharacterType;

    fn index(&self, index: usize) -> &Self::Output {
        self.as_genome_subsequence().index(index)
    }
}

impl<AlphabetType: Alphabet> Index<Range<usize>> for VectorSubGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        VectorSubGenome::ref_cast(&self.slice[index.start..index.end])
    }
}

impl<AlphabetType: Alphabet> Index<RangeFrom<usize>> for VectorSubGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeFrom<usize>) -> &Self::Output {
        self.index(index.start..self.len())
    }
}

impl<AlphabetType: Alphabet> Index<RangeTo<usize>> for VectorSubGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeTo<usize>) -> &Self::Output {
        self.index(0..index.end)
    }
}

impl<AlphabetType: Alphabet> Index<RangeFull> for VectorSubGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, _index: RangeFull) -> &Self::Output {
        self.index(0..self.len())
    }
}

impl<AlphabetType: Alphabet> Index<RangeInclusive<usize>> for VectorSubGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeInclusive<usize>) -> &Self::Output {
        VectorSubGenome::ref_cast(&self.slice[*index.start()..=*index.end()])
    }
}

impl<AlphabetType: Alphabet> Index<RangeToInclusive<usize>> for VectorSubGenome<AlphabetType> {
    type Output = VectorSubGenome<AlphabetType>;

    fn index(&self, index: RangeToInclusive<usize>) -> &Self::Output {
        self.index(0..=index.end)
    }
}

impl<AlphabetType: Alphabet> Index<usize> for VectorSubGenome<AlphabetType> {
    type Output = AlphabetType::CharacterType;

    fn index(&self, index: usize) -> &Self::Output {
        self.slice.index(index)
    }
}

impl<AlphabetType: Alphabet> FromIterator<AlphabetType::CharacterType>
    for VectorGenome<AlphabetType>
{
    fn from_iter<T: IntoIterator<Item = AlphabetType::CharacterType>>(iter: T) -> Self {
        let mut result = Self::default();
        result.extend(iter);
        result
    }
}

impl<AlphabetType: Alphabet> Extend<AlphabetType::CharacterType> for VectorGenome<AlphabetType> {
    fn extend<T: IntoIterator<Item = AlphabetType::CharacterType>>(&mut self, iter: T) {
        self.vector.extend(iter)
    }
}

impl<AlphabetType: Alphabet> IntoIterator for VectorGenome<AlphabetType> {
    type Item = AlphabetType::CharacterType;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.vector.into_iter()
    }
}

impl<AlphabetType: Alphabet> Borrow<VectorSubGenome<AlphabetType>> for VectorGenome<AlphabetType> {
    fn borrow(&self) -> &VectorSubGenome<AlphabetType> {
        self.as_genome_subsequence()
    }
}

impl<AlphabetType: Alphabet> ToOwned for VectorSubGenome<AlphabetType> {
    type Owned = VectorGenome<AlphabetType>;

    fn to_owned(&self) -> Self::Owned {
        self.iter().cloned().collect()
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    SequenceMut<'a, AlphabetType::CharacterType, VectorSubGenome<AlphabetType>>
    for VectorGenome<AlphabetType>
{
    type IteratorMut = std::slice::IterMut<'a, AlphabetType::CharacterType>;

    fn iter_mut(&'a mut self) -> Self::IteratorMut {
        self.vector.iter_mut()
    }
}

impl<'a, AlphabetType: Alphabet + 'a> IndexMut<usize> for VectorGenome<AlphabetType> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.vector.index_mut(index)
    }
}

impl<'a, AlphabetType: Alphabet + 'a> IndexMut<Range<usize>> for VectorGenome<AlphabetType> {
    fn index_mut(&mut self, index: Range<usize>) -> &mut Self::Output {
        VectorSubGenome::ref_cast_mut(&mut self.vector[index])
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    SequenceMut<'a, AlphabetType::CharacterType, VectorSubGenome<AlphabetType>>
    for VectorSubGenome<AlphabetType>
{
    type IteratorMut = std::slice::IterMut<'a, AlphabetType::CharacterType>;

    fn iter_mut(&'a mut self) -> Self::IteratorMut {
        self.slice.iter_mut()
    }
}

impl<'a, AlphabetType: Alphabet + 'a> IndexMut<usize> for VectorSubGenome<AlphabetType> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.slice.index_mut(index)
    }
}

impl<'a, AlphabetType: Alphabet + 'a> IndexMut<Range<usize>> for VectorSubGenome<AlphabetType> {
    fn index_mut(&mut self, index: Range<usize>) -> &mut Self::Output {
        VectorSubGenome::ref_cast_mut(&mut self.slice[index])
    }
}

impl<AlphabetType: Alphabet> Default for VectorGenome<AlphabetType> {
    fn default() -> Self {
        Self {
            vector: Default::default(),
        }
    }
}

impl<'a, AlphabetType: Alphabet + 'a>
    OwnedSequence<'a, AlphabetType::CharacterType, VectorSubGenome<AlphabetType>>
    for VectorGenome<AlphabetType>
{
}

#[cfg(test)]
mod tests {
    use crate::implementation::vec_sequence::VectorGenome;
    use crate::interface::alphabet::dna_alphabet::DnaAlphabet;
    use crate::interface::sequence::{GenomeSequence, OwnedGenomeSequence};

    #[test]
    fn test_reverse_complement() {
        let genome = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();
        let reverse_complement = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ACCGAAT").unwrap();
        debug_assert_eq!(genome.clone_as_reverse_complement(), reverse_complement);
        debug_assert_eq!(genome, reverse_complement.clone_as_reverse_complement());
    }

    #[test]
    fn test_display() {
        let genome = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ATTCGGT").unwrap();
        let display_string = genome.as_string();
        let expected_string = "ATTCGGT";
        debug_assert_eq!(display_string, expected_string);
    }
}
