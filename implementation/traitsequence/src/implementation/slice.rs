use crate::interface::{CloneableSequence, Sequence, SequenceMut};

impl<'a, Item: 'a> Sequence<'a, Item, [Item]> for [Item] {
    type Iterator = std::slice::Iter<'a, Item>;

    fn iter(&'a self) -> Self::Iterator {
        <[Item]>::iter(self)
    }

    fn len(&self) -> usize {
        <[Item]>::len(self)
    }
}

impl<'a, Item: 'a> SequenceMut<'a, Item, [Item]> for [Item] {
    type IteratorMut = std::slice::IterMut<'a, Item>;
    fn iter_mut(&'a mut self) -> Self::IteratorMut {
        <[Item]>::iter_mut(self)
    }
}

impl<'a, Item: 'a + Clone> CloneableSequence<'a, Item, [Item]> for [Item] {}

#[cfg(test)]
mod tests {
    use crate::interface::Sequence;

    #[test]
    fn test_len() {
        let array = [0, 1, 2];
        let slice = &array[0..2];
        // Making sure that the fully qualified syntax in the trait implementation works as I think and does not create endless recursion.
        debug_assert_eq!(2, Sequence::len(slice));
    }
}
