use std::iter::FromIterator;

/// A genome string.
/// It should be lexically ordered.
pub trait Genome: for<'a> FromIterator<&'a u8> + FromIterator<u8> + Eq + Clone + Ord
where
    for<'a> &'a Self: IntoIterator<Item = u8>,
{
    /**
     * Returns the reverse complement of this genome.
     *
     * Panics if this genome is [not valid](is_valid).
     */
    fn reverse_complement(&self) -> Self;

    /**
     * Returns true if this genome is valid, i.e. it contains no invalid characters.
     *
     * Valid characters are defined by [is_valid_ascii_genome_character()](is_valid_ascii_genome_character)
     */
    fn is_valid(&self) -> bool {
        self.into_iter().all(is_valid_ascii_genome_character)
    }

    /// Returns the amount of bases in this genome string.
    fn len(&self) -> usize;

    /// Returns true if this genome string is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Copies this genome string into a `Vec`.
    fn into_vec(&self) -> Vec<u8> {
        self.into_iter().collect()
    }
}

/**
 * Returns the complement of the given genome char.
 * Returns `None` if the given char [is invalid](is_valid_ascii_genome_character).
 */
pub fn ascii_complement(char: u8) -> Option<u8> {
    match char {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'G' => Some(b'C'),
        b'C' => Some(b'G'),
        _ => None,
    }
}

/**
 * Returns true if the given ascii character represents a valid genome character.
 * Valid genome characters are `A`, `T`, `G` and `C`.
 */
pub fn is_valid_ascii_genome_character(char: u8) -> bool {
    match char {
        b'A' | b'T' | b'G' | b'C' => true,
        _ => false,
    }
}
