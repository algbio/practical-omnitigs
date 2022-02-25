pub mod vec_sequence;

pub mod bit_vec_sequence;
pub mod bit_vec_sequence_store;
pub mod vec_sequence_store;

/// The default genome type that achieves a good balance between speed and size.
pub type DefaultGenome<AlphabetType> = bit_vec_sequence::BitVectorGenome<AlphabetType>;
/// The default genome subsequence type that achieves a good balance between speed and size.
pub type DefaultSubGenome<AlphabetType> = bit_vec_sequence::BitVectorSubGenome<AlphabetType>;
/// The default genome sequence store type that achieves a good balance between speed and size.
pub type DefaultSequenceStore<AlphabetType> =
    bit_vec_sequence_store::BitVectorSequenceStore<AlphabetType>;
/// The handle type of the default genome sequence store type.
pub type DefaultSequenceStoreHandle<AlphabetType> = <bit_vec_sequence_store::BitVectorSequenceStore<
    AlphabetType,
> as crate::interface::sequence_store::SequenceStore<AlphabetType>>::Handle;
