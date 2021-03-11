use crate::interface::{Sequence, SequenceMut};

impl<'a, Item: 'a> Sequence<'a, Item> for &'a [Item] {
    type Iterator = std::slice::Iter<'a, Item>;

    fn iter(self) -> Self::Iterator {
        self[..].iter()
    }

    fn len(self) -> usize {
        Vec::len(self)
    }
}

impl<'a, Item: 'a> SequenceMut<'a, Item> for &'a mut Vec<Item> {
    type IteratorMut = std::slice::IterMut<'a, Item>;

    fn iter_mut(mut self) -> Self::IteratorMut {
        self[..].iter_mut()
    }
}