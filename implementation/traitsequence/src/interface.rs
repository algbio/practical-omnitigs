use std::fmt::{Debug, Write};
use std::iter::FromIterator;
use std::ops::{Index, IndexMut, Range};

/// A type behaving like a sequence over the type `Item`.
pub trait Sequence<'a, Item: 'a, Subsequence: Sequence<'a, Item, Subsequence> + ?Sized>:
    Index<usize, Output = Item> + Index<Range<usize>, Output = Subsequence>
{
    /// The iterator type of the sequence.
    type Iterator: DoubleEndedIterator<Item = &'a Item>;

    /// Returns a prefix with length `len` of this sequence.
    /// Panics if `len >= self.len()`.
    fn prefix(&'a self, len: usize) -> &Subsequence {
        debug_assert!(len < self.len());
        &self[0..len]
    }

    /// Returns a suffix with length `len` of this sequence.
    /// Panics if `len >= self.len()`.
    fn suffix(&'a self, len: usize) -> &Subsequence {
        debug_assert!(len < self.len());
        &self[self.len() - len..self.len()]
    }

    /// Returns an iterator over the sequence.
    fn iter(&'a self) -> Self::Iterator;

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

    /// Converts the sequence to a string using the debug formatting of the items.
    ///
    /// ```rust
    /// use traitsequence::interface::Sequence;
    ///
    /// let sequence = [0, 2, 1];
    /// debug_assert_eq!(sequence.to_debug_string(), "[0, 2, 1]".to_string());
    ///
    /// let sequence = ["a", "c", "b"];
    /// debug_assert_eq!(sequence.to_debug_string(), "[\"a\", \"c\", \"b\"]".to_string());
    /// ```
    fn to_debug_string(&'a self) -> String
    where
        Item: Debug,
    {
        let mut result = String::new();
        write!(result, "[").unwrap();
        let mut once = true;
        for item in self.iter() {
            if once {
                once = false;
            } else {
                write!(result, ", ").unwrap();
            }
            write!(result, "{:?}", item).unwrap();
        }
        write!(result, "]").unwrap();
        result
    }
}

/// A type behaving like a mutable sequence over the type `Item`.
/// That is, its items can be mutated, but the sequence it self can not.
/// For a sequence where items can be appended, rearranged etc. see [EditableSequence].
pub trait SequenceMut<'a, Item: 'a, Subsequence: SequenceMut<'a, Item, Subsequence> + ?Sized>:
    Sequence<'a, Item, Subsequence>
    + IndexMut<usize, Output = Item>
    + IndexMut<Range<usize>, Output = Subsequence>
{
    /// The mutable iterator type of the sequence.
    type IteratorMut: Iterator<Item = &'a mut Item>;

    /// Returns a mutable iterator over the sequence.
    fn iter_mut(&'a mut self) -> Self::IteratorMut;
}

/// A type behaving like an owned sequence over the type `Item`.
/// Currently this only means the sequence is `Sized`.
pub trait OwnedSequence<'a, Item: 'a, Subsequence: Sequence<'a, Item, Subsequence> + ?Sized>:
    Sequence<'a, Item, Subsequence> + Sized
{
}

/// A type behaving like an cloneable sequence over the type `Item`.
/// Currently this only means the sequence is `ToOwned`.
pub trait CloneableSequence<
    'a,
    Item: 'a + Clone,
    Subsequence: CloneableSequence<'a, Item, Subsequence> + ?Sized,
>: ToOwned
{
}

/// A type behaving like a sequence over the type `Item` that can be edited.
/// This sequences items can not necessarily be mutated themselves, but they can be rearranged or new items can be appended etc.
/// For a sequence where the items themselves can be mutated, see [SequenceMut].
pub trait EditableSequence<'a, Item: 'a, Subsequence: Sequence<'a, Item, Subsequence> + ?Sized>:
    Sequence<'a, Item, Subsequence> + Extend<Item> + IntoIterator<Item = Item> + FromIterator<Item>
{
    /// Extend this sequence from a sequence of compatible items.
    fn extend_into<
        ExtensionItem: Into<Item>,
        ExtensionSource: IntoIterator<Item = ExtensionItem>,
    >(
        &mut self,
        extension: ExtensionSource,
    ) {
        self.extend(extension.into_iter().map(Into::into));
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
        debug_assert_eq!(merged, vec![0, 1, 2, 3, 4, 5, 6, 7, 8]);
    }
}
