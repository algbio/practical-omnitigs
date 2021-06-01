pub mod ascii_vec_sequence;
pub mod two_bit_vec_sequence;

pub mod two_bit_vec_sequence_store;
pub mod vec_sequence_store;

/// The default genome type that achieves a good balance between speed and size.
pub type DefaultGenome = two_bit_vec_sequence::TwoBitVectorGenome;
/// The default genome subsequence type that achieves a good balance between speed and size.
pub type DefaultSubGenome = two_bit_vec_sequence::TwoBitVectorSubGenome;
/// The default genome sequence store type that achieves a good balance between speed and size.
pub type DefaultSequenceStore = two_bit_vec_sequence_store::TwoBitVectorSequenceStore;
/// The handle type of the default genome sequence store type.
pub type DefaultSequenceStoreHandle = <two_bit_vec_sequence_store::TwoBitVectorSequenceStore as crate::interface::sequence_store::SequenceStore>::Handle;
