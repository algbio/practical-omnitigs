use crate::ascii_complement;
use crate::genome::Genome;
use std::iter::{Cloned, FromIterator};

pub struct VectorGenome(Vec<u8>);

impl Genome for VectorGenome {
    fn reverse_complement(&self) -> Self {
        self.0
            .iter()
            .rev()
            .cloned()
            .map(|c| {
                ascii_complement(c)
                    .expect("Reverse complement can only be computed from valid genome strings")
            })
            .collect()
    }
}

impl std::fmt::Display for VectorGenome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for c in self {
            write!(f, "{}", c as char)?;
        }

        Ok(())
    }
}

impl<'a> IntoIterator for &'a VectorGenome {
    type Item = u8;
    type IntoIter = Cloned<std::slice::Iter<'a, u8>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter().cloned()
    }
}

impl FromIterator<u8> for VectorGenome {
    fn from_iter<T: IntoIterator<Item = u8>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}
