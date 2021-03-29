//! A simple representation of a genome as `Vec<u8>`.

use crate::interface::ascii_complement;
use crate::interface::{ExtendableGenome, Genome};
use std::iter::{Cloned, FromIterator};
use traitsequence::interface::Sequence;

/// A simple representation of a genome as `Vec<u8>`.
/// This is not very efficient, but was quick to implement.
#[derive(Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Hash, Default)]
pub struct VectorGenome(Vec<u8>);

impl VectorGenome {
    /// Creates a new VectorGenome from the given characters.
    pub fn new(genome: Vec<u8>) -> Self {
        Self(genome)
    }
}

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

impl ExtendableGenome for VectorGenome {
    fn extend<ExtensionSource: IntoIterator<Item = u8>>(&mut self, extension: ExtensionSource) {
        self.0.extend(extension)
    }
}

impl<'a> Sequence<'a, u8> for VectorGenome {
    type Iterator = std::slice::Iter<'a, u8>;
    type IteratorMut = std::slice::IterMut<'a, u8>;

    fn iter(&'a self) -> Self::Iterator {
        self.0.iter()
    }

    fn iter_mut(&'a mut self) -> Self::IteratorMut {
        self.0.iter_mut()
    }

    fn len(&self) -> usize {
        self.0.len()
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

impl<'a> FromIterator<&'a u8> for VectorGenome {
    fn from_iter<T: IntoIterator<Item = &'a u8>>(iter: T) -> Self {
        Self(iter.into_iter().cloned().collect())
    }
}

impl FromIterator<u8> for VectorGenome {
    fn from_iter<T: IntoIterator<Item = u8>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<IndexType> std::ops::Index<IndexType> for VectorGenome
where
    Vec<u8>: std::ops::Index<IndexType>,
{
    type Output = <Vec<u8> as std::ops::Index<IndexType>>::Output;

    fn index(&self, index: IndexType) -> &Self::Output {
        self.0.index(index)
    }
}

#[cfg(test)]
mod tests {
    use crate::implementation::vector_genome_impl::VectorGenome;
    use crate::interface::Genome;
    use std::iter::FromIterator;

    #[test]
    fn test_reverse_complement() {
        let genome = VectorGenome::from_iter(b"ATTCGGT");
        let reverse_complement = VectorGenome::from_iter(b"ACCGAAT");
        assert_eq!(genome.reverse_complement(), reverse_complement);
        assert_eq!(genome, reverse_complement.reverse_complement());
    }

    #[test]
    fn test_display() {
        let genome = VectorGenome::from_iter(b"ATTCGGT");
        let display_string = genome.to_string();
        let expected_string = "ATTCGGT";
        assert_eq!(display_string, expected_string);
    }
}
