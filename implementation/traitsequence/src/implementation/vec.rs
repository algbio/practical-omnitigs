use crate::interface::{CloneableSequence, EditableSequence, OwnedSequence, Sequence, SequenceMut};

impl<'a, Item: 'a> Sequence<'a, Item, [Item]> for Vec<Item> {
    type Iterator = std::slice::Iter<'a, Item>;
    fn iter(&'a self) -> Self::Iterator {
        self[..].iter()
    }

    fn len(&self) -> usize {
        Vec::len(self)
    }
}

impl<'a, Item: 'a> SequenceMut<'a, Item, [Item]> for Vec<Item> {
    type IteratorMut = std::slice::IterMut<'a, Item>;

    fn iter_mut(&'a mut self) -> Self::IteratorMut {
        self[..].iter_mut()
    }
}

impl<'a, Item: 'a> OwnedSequence<'a, Item, [Item]> for Vec<Item> {}

impl<'a, Item: 'a + Clone> CloneableSequence<'a, Item, [Item]> for Vec<Item> {}

impl<'a, Item: 'a> EditableSequence<'a, Item, [Item]> for Vec<Item> {}
