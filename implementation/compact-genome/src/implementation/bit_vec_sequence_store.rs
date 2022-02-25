//! A vector sequence store where each character is encoded as two bits.

use crate::implementation::bit_vec_sequence::{
    alphabet_character_bit_width, BitVectorGenome, BitVectorSubGenome,
};
use crate::interface::alphabet::{Alphabet, AlphabetCharacter, AlphabetError};
use crate::interface::sequence::GenomeSequence;
use crate::interface::sequence_store::{
    HandleWithLength, HandleWithSubsequence, InverseMappingSequenceStore, SequenceStore,
};
use bitvec::order::Lsb0;
use bitvec::view::BitView;
use std::marker::PhantomData;
use traitsequence::interface::Sequence;

/// A bitvector based sequence store.
#[derive(Default, Clone, Eq, PartialEq, Debug)]
pub struct BitVectorSequenceStore<AlphabetType: Alphabet> {
    sequence: BitVectorGenome<AlphabetType>,
}

/// A handle of a sequence in an [BitVectorSequenceStore].
#[derive(Default, Debug, Clone, Copy, Eq, PartialEq)]
pub struct BitVectorSequenceStoreHandle<AlphabetType: Alphabet> {
    offset: usize,
    len: usize,
    phantom_data: PhantomData<AlphabetType>,
}

impl<AlphabetType: Alphabet> BitVectorSequenceStore<AlphabetType> {
    /// Creates a new instance.
    pub fn new() -> Self {
        Self {
            sequence: Default::default(),
        }
    }

    /// Returns the number of bytes consumed by the characters stored in this sequence store.
    pub fn size_in_memory(&self) -> usize {
        (self.sequence.len() - 1) / 4 + 1 // Rounding up integer division
    }
}

impl<AlphabetType: Alphabet + 'static> SequenceStore<AlphabetType>
    for BitVectorSequenceStore<AlphabetType>
{
    type Handle = BitVectorSequenceStoreHandle<AlphabetType>;
    type SequenceRef = BitVectorSubGenome<AlphabetType>;

    fn add<
        Sequence: for<'a> GenomeSequence<'a, AlphabetType, Subsequence> + ?Sized,
        Subsequence: for<'a> GenomeSequence<'a, AlphabetType, Subsequence> + ?Sized,
    >(
        &mut self,
        s: &Sequence,
    ) -> Self::Handle {
        let offset = self.sequence.len();
        let len = s.len();
        self.sequence.extend(s.iter().cloned());
        Self::Handle {
            offset,
            len,
            phantom_data: Default::default(),
        }
    }

    fn add_from_iter_u8<IteratorType: IntoIterator<Item = u8>>(
        &mut self,
        iter: IteratorType,
    ) -> Result<Self::Handle, AlphabetError> {
        let offset = self.sequence.len();
        let iter = iter.into_iter();
        let (size, _) = iter.size_hint();
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        self.sequence.bits.reserve(size * bit_width);
        for item in iter {
            match AlphabetType::ascii_to_character(item) {
                Ok(character) => self
                    .sequence
                    .bits
                    .extend_from_bitslice(&character.index().view_bits::<Lsb0>()[0..bit_width]),

                Err(error) => {
                    self.sequence.bits.resize(offset * bit_width, false);
                    return Err(error);
                }
            }
        }

        let len = self.sequence.len() - offset;
        Ok(Self::Handle {
            offset,
            len,
            phantom_data: Default::default(),
        })
    }

    fn get(&self, handle: &Self::Handle) -> &Self::SequenceRef {
        &self.sequence[handle.offset..handle.offset + handle.len]
    }
}

impl<AlphabetType: Alphabet + 'static> InverseMappingSequenceStore<AlphabetType>
    for BitVectorSequenceStore<AlphabetType>
{
    fn map_sequence_ref_to_handle(&self, sequence_ref: &Self::SequenceRef) -> Self::Handle {
        let raw_offset = self.sequence.bits.offset_from(&sequence_ref.bits[..]);
        debug_assert!(raw_offset >= 0);
        let bit_width = alphabet_character_bit_width(AlphabetType::SIZE);
        let offset = raw_offset as usize / bit_width;

        Self::Handle {
            offset,
            len: sequence_ref.len(),
            phantom_data: Default::default(),
        }
    }
}

impl<AlphabetType: Alphabet> HandleWithLength for BitVectorSequenceStoreHandle<AlphabetType> {
    fn len(&self) -> usize {
        self.len
    }
}

impl<AlphabetType: Alphabet> HandleWithSubsequence<core::ops::Range<usize>>
    for BitVectorSequenceStoreHandle<AlphabetType>
{
    fn subsequence_handle(&self, range: core::ops::Range<usize>) -> Self {
        let result = Self {
            offset: self.offset + range.start,
            len: range.end - range.start,
            phantom_data: self.phantom_data,
        };
        debug_assert!(self.offset + self.len >= result.offset + result.len);
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::implementation::bit_vec_sequence_store::BitVectorSequenceStore;
    use crate::implementation::vec_sequence::VectorGenome;
    use crate::interface::alphabet::dna_alphabet::DnaAlphabet;
    use crate::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
    use crate::interface::sequence_store::{InverseMappingSequenceStore, SequenceStore};

    #[test]
    fn test_inverse_mapping() {
        let mut sequence_store = BitVectorSequenceStore::<DnaAlphabet>::new();
        let handle1 = sequence_store.add_from_slice_u8(b"ACGTTG").unwrap();
        let handle2 = sequence_store.add_from_slice_u8(b"CGACTG").unwrap();
        let reference1 = sequence_store.get(&handle1);
        let reference2 = sequence_store.get(&handle2);
        debug_assert_eq!(
            reference1.convert::<VectorGenome<_>, _>(),
            VectorGenome::from_slice_u8(b"ACGTTG").unwrap()
        );
        debug_assert_eq!(
            reference2.convert::<VectorGenome<_>, _>(),
            VectorGenome::from_slice_u8(b"CGACTG").unwrap()
        );
        debug_assert_eq!(
            sequence_store.map_sequence_ref_to_handle(reference1),
            handle1
        );
        debug_assert_eq!(
            sequence_store.map_sequence_ref_to_handle(reference2),
            handle2
        );
    }
}
