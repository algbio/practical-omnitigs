/// A type behaving like a sequence over the type `Item` with a subsequence type `Subsequence`.
pub trait Sequence<'a, Item: 'a>: std::ops::Index<usize, Output = Item> {
    /// The iterator type of the sequence.
    type Iterator: Iterator<Item = &'a Item>;
    /// The mutable iterator type of the sequence.
    type IteratorMut: Iterator<Item = &'a mut Item>;

    /// Returns an iterator over the sequence.
    fn iter(&'a self) -> Self::Iterator;

    /// Returns a mutable iterator over the sequence.
    fn iter_mut(&'a mut self) -> Self::IteratorMut;

    /// Returns the length of the sequence.
    fn len(&self) -> usize;

    /// Returns true if the sequence is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the first item of the sequence.
    fn first(&'a self) -> Option<&Item> {
        self.iter().next()
    }

    /// Returns the last item of the sequence.
    fn last(&'a self) -> Option<&Item> {
        self.iter().last()
    }

    /// Returns true if this is a proper subsequence of the given sequence.
    /// Proper means that the sequences are not equal.
    fn is_proper_subsequence_of(&'a self, other: &Self) -> bool
    where
        Item: Eq,
    {
        if self.len() >= other.len() {
            return false;
        }

        for start_index in 0..=other.len() - self.len() {
            let mut found_subsequence = true;
            for index in 0..self.len() {
                if self[index] != other[start_index + index] {
                    found_subsequence = false;
                    break;
                }
            }
            if found_subsequence {
                return true;
            }
        }

        false
    }
}
