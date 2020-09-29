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
}
