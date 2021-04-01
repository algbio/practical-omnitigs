//! A simple represenation of a genome as `[u8]`

use crate::interface::{GenomeSequence, GenomeSequenceMut};

impl<'a> GenomeSequence<'a, [u8]> for [u8] {}

impl<'a> GenomeSequenceMut<'a, [u8]> for [u8] {}
