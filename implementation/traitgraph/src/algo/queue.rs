/// A queue that supports both popping and pushing at front and back.
pub trait BidirectedQueue<T>: Default {
    /// Insert an element at the front of the queue.
    fn push_front(&mut self, t: T);
    /// Insert an element at the back of the queue.
    fn push_back(&mut self, t: T);
    /// Remove and return an element from the front of the queue.
    fn pop_front(&mut self) -> Option<T>;
    /// Remove and return an element from the back of the queue.
    fn pop_back(&mut self) -> Option<T>;
    /// Remove all elements from the queue without returning them.
    fn clear(&mut self);
    /// Return the amount of elements currently in the queue.
    fn len(&self) -> usize;
    /// Returns true if the queue contains no elements.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<T> BidirectedQueue<T> for std::collections::LinkedList<T> {
    fn push_front(&mut self, t: T) {
        std::collections::LinkedList::<T>::push_front(self, t)
    }

    fn push_back(&mut self, t: T) {
        std::collections::LinkedList::<T>::push_back(self, t)
    }

    fn pop_front(&mut self) -> Option<T> {
        std::collections::LinkedList::<T>::pop_front(self)
    }

    fn pop_back(&mut self) -> Option<T> {
        std::collections::LinkedList::<T>::pop_back(self)
    }

    fn clear(&mut self) {
        self.clear();
    }

    fn len(&self) -> usize {
        std::collections::LinkedList::<T>::len(self)
    }
}

impl<T> BidirectedQueue<T> for std::collections::VecDeque<T> {
    fn push_front(&mut self, t: T) {
        std::collections::VecDeque::<T>::push_front(self, t)
    }

    fn push_back(&mut self, t: T) {
        std::collections::VecDeque::<T>::push_back(self, t)
    }

    fn pop_front(&mut self) -> Option<T> {
        std::collections::VecDeque::<T>::pop_front(self)
    }

    fn pop_back(&mut self) -> Option<T> {
        std::collections::VecDeque::<T>::pop_back(self)
    }

    fn clear(&mut self) {
        self.clear();
    }

    fn len(&self) -> usize {
        std::collections::VecDeque::<T>::len(self)
    }
}
