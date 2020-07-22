use crate::{ImmutableGraphContainer, NavigableGraph, NodeIndex};
use num_traits::{NumCast, PrimInt};
use std::collections::LinkedList;

pub struct Bfs<IndexType> {
    queue: LinkedList<NodeIndex<IndexType>>,
    order: Vec<IndexType>,
    current_order: IndexType,
}

impl<IndexType: PrimInt> Bfs<IndexType> {
    pub fn new<
        NodeData,
        EdgeData,
        Graph: ImmutableGraphContainer<NodeData, EdgeData, IndexType>,
    >(
        graph: &Graph,
        start: NodeIndex<IndexType>,
    ) -> Self {
        let mut queue = LinkedList::new();
        queue.push_back(start);
        let mut order = vec![IndexType::max_value(); graph.node_count()];
        order[start] = IndexType::zero();
        Self {
            queue,
            order,
            current_order: IndexType::one(),
        }
    }

    pub fn next<
        'a,
        NodeData,
        EdgeData,
        Graph: ImmutableGraphContainer<NodeData, EdgeData, IndexType>
            + NavigableGraph<'a, NodeData, EdgeData, IndexType>,
    >(
        &mut self,
        graph: &'a Graph,
    ) -> Option<NodeIndex<IndexType>> {
        if let Some(first) = self.queue.pop_front() {
            for neighbor in graph.out_neighbors(first).unwrap() {
                let order_entry =
                    &mut self.order[<usize as NumCast>::from(neighbor.node_id).unwrap()];
                if *order_entry == IndexType::max_value() {
                    *order_entry = self.current_order;
                    self.current_order = self.current_order + IndexType::one();
                    self.queue.push_back(neighbor.node_id);
                }
            }

            Some(first)
        } else {
            None
        }
    }

    pub fn order_of(&self, node: NodeIndex<IndexType>) -> Option<IndexType> {
        let order = self.order[<usize as NumCast>::from(node).unwrap()];
        if order == IndexType::max_value() {
            None
        } else {
            Some(order)
        }
    }
}
