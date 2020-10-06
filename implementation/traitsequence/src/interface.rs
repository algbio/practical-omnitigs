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

    /// Returns true if this sequence contains the given item.
    fn contains(&'a self, item: &Item) -> bool
    where
        Item: Eq,
    {
        self.iter().any(|i| item == i)
    }

    /// Returns an iterator over this sequence merged before the given other sequence under the assumption that the sequences can be merged this way.
    /// A merge is possible if a non-empty suffix of this sequence equals a non-empty prefix of the other sequence.
    ///
    /// The method panics if this sequence does not contain the first item of the other sequence or the other sequence is empty.
    /// The method does not fail if the sequences are not mergeable for other reasons.
    fn forward_merge_iter_assume_mergeable(
        &'a self,
        suffix: &'a Self,
    ) -> std::iter::Chain<Self::Iterator, std::iter::Skip<Self::Iterator>>
    where
        Item: Eq,
    {
        let first_item_index = self
            .iter()
            .enumerate()
            .filter(|(_, i)| *i == suffix.first().expect("The given sequence is empty."))
            .map(|(i, _)| i)
            .next()
            .expect("This sequence does not contain the first item of the given sequence.");
        self.iter()
            .chain(suffix.iter().skip(self.len() - first_item_index))
    }

    /// Returns an iterator over this sequence merged after the given other sequence under the assumption that the sequences can be merged this way.
    /// A merge is possible if a non-empty prefix of this sequence equals a non-empty suffix of the other sequence.
    ///
    /// The method panics if the other sequence does not contain the first item of this sequence or this sequence is empty.
    /// The method does not fail if the sequences are not mergeable for other reasons.
    fn backward_merge_iter_assume_mergeable(
        &'a self,
        suffix: &'a Self,
    ) -> std::iter::Chain<Self::Iterator, std::iter::Skip<Self::Iterator>>
    where
        Item: Eq,
    {
        suffix.forward_merge_iter_assume_mergeable(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::interface::Sequence;

    #[test]
    fn test_merge_sequences_simple() {
        let s1 = vec![0, 1, 2, 3, 4, 5];
        let s2 = vec![3, 4, 5, 6, 7, 8];
        let merged: Vec<_> = s1
            .forward_merge_iter_assume_mergeable(&s2)
            .copied()
            .collect();
        assert_eq!(merged, vec![0, 1, 2, 3, 4, 5, 6, 7, 8]);
    }
}
