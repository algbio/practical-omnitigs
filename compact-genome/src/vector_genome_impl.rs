use crate::genome::Genome;
use std::fmt::Formatter;

pub struct VectorGenome(Vec<u8>);

impl Genome<'_> for VectorGenome {
    fn reverse_complement(&self) -> Self {
        unimplemented!()
    }
}

impl std::fmt::Display for VectorGenome {
    fn fmt(&self, _f: &mut Formatter<'_>) -> std::fmt::Result {
        unimplemented!()
    }
}

impl From<&[u8]> for VectorGenome {
    fn from(slice: &[u8]) -> Self {
        Self (slice.into())
    }
}

impl From<VectorGenome> for Vec<u8> {
    fn from(genome: VectorGenome) -> Self {
        genome.0
    }
}