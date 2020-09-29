use crate::interface::Sequence;

impl<'a, Item: 'a> Sequence<'a, Item> for Vec<Item> {
    type Iterator = std::slice::Iter<'a, Item>;
    type IteratorMut = std::slice::IterMut<'a, Item>;

    fn iter(&'a self) -> Self::Iterator {
        self[..].iter()
    }

    fn iter_mut(&'a mut self) -> Self::IteratorMut {
        self[..].iter_mut()
    }

    fn len(&self) -> usize {
        Vec::len(self)
    }
}
