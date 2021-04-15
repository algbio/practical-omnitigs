pub mod ascii_vec_sequence;
pub mod slice_sequence;
pub mod two_bit_vec_sequence;

pub mod vec_sequence_store;

/// The default genome type that achieves a good balance between speed and size.
pub type DefaultGenome = two_bit_vec_sequence::TwoBitVectorGenome;
/// The default genome subsequence type that achieves a good balance between speed and size.
pub type DefaultSubGenome = two_bit_vec_sequence::TwoBitVectorSubGenome;
