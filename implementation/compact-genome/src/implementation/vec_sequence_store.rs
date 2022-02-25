//! An ASCII vector based sequence store.

use crate::implementation::vec_sequence::VectorSubGenome;
use crate::interface::alphabet::{Alphabet, AlphabetCharacter, AlphabetError};
use crate::interface::sequence::GenomeSequence;
use crate::interface::sequence_store::{
    HandleWithLength, HandleWithSubsequence, InverseMappingSequenceStore, SequenceStore,
};
use ref_cast::RefCast;
use std::marker::PhantomData;
use traitsequence::interface::Sequence;

/// An plain vector based sequence store.
#[derive(Default, Clone, Eq, PartialEq, Debug)]
pub struct VectorSequenceStore<AlphabetType: Alphabet> {
    sequence: Vec<AlphabetType::CharacterType>,
}

/// A handle of a sequence in an [VectorSequenceStore].
#[derive(Default, Clone, Copy, Debug, Eq, PartialEq)]
pub struct VectorSequenceStoreHandle<AlphabetType: Alphabet> {
    offset: usize,
    len: usize,
    phantom_data: PhantomData<AlphabetType>,
}

impl<AlphabetType: Alphabet> VectorSequenceStore<AlphabetType> {
    /// Creates a new instance.
    pub fn new() -> Self {
        Self {
            sequence: Vec::new(),
        }
    }
}

impl<AlphabetType: Alphabet + 'static> SequenceStore<AlphabetType>
    for VectorSequenceStore<AlphabetType>
{
    type Handle = VectorSequenceStoreHandle<AlphabetType>;
    type SequenceRef = VectorSubGenome<AlphabetType>;

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
        self.sequence.reserve(size);
        for item in iter {
            match AlphabetType::ascii_to_character(item) {
                Ok(character) => self.sequence.push(character),
                Err(error) => {
                    self.sequence
                        .resize(offset, AlphabetType::CharacterType::from_index(0).unwrap());
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
        VectorSubGenome::ref_cast(&self.sequence[handle.offset..handle.offset + handle.len])
    }
}

impl<AlphabetType: Alphabet + 'static> InverseMappingSequenceStore<AlphabetType>
    for VectorSequenceStore<AlphabetType>
{
    fn map_sequence_ref_to_handle(&self, sequence_ref: &Self::SequenceRef) -> Self::Handle {
        Self::Handle {
            offset: (sequence_ref.slice.as_ptr() as usize) - (self.sequence.as_ptr() as usize),
            len: sequence_ref.len(),
            phantom_data: Default::default(),
        }
    }
}

impl<AlphabetType: Alphabet> HandleWithLength for VectorSequenceStoreHandle<AlphabetType> {
    fn len(&self) -> usize {
        self.len
    }
}

impl<AlphabetType: Alphabet> HandleWithSubsequence<core::ops::Range<usize>>
    for VectorSequenceStoreHandle<AlphabetType>
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
