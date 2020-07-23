use crate::algo::queue::BidirectedQueue;
use crate::{NavigableGraph, NodeIndex, StaticGraph};
use num_traits::{NumCast, PrimInt};
use std::collections::LinkedList;
use std::marker::PhantomData;

/// A normal forward BFS in a directed graph.
pub type ForwardBfs<NodeData, EdgeData, IndexType, Graph> = Traversal<
    NodeData,
    EdgeData,
    IndexType,
    Graph,
    ForwardNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<NodeIndex<IndexType>>,
>;
/// A normal backward BFS in a directed graph.
pub type BackwardBfs<NodeData, EdgeData, IndexType, Graph> = Traversal<
    NodeData,
    EdgeData,
    IndexType,
    Graph,
    BackwardNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<NodeIndex<IndexType>>,
>;
/// A BFS that treats each directed edge as an undirected edge, i.e. that traverses edge both in forward and backward direction.
pub type UndirectedBfs<NodeData, EdgeData, IndexType, Graph> = Traversal<
    NodeData,
    EdgeData,
    IndexType,
    Graph,
    UndirectedNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<NodeIndex<IndexType>>,
>;

/// A normal forward DFS in a directed graph.
pub type ForwardDfs<NodeData, EdgeData, IndexType, Graph> = Traversal<
    NodeData,
    EdgeData,
    IndexType,
    Graph,
    ForwardNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<NodeIndex<IndexType>>,
>;
/// A normal backward DFS in a directed graph.
pub type BackwardDfs<NodeData, EdgeData, IndexType, Graph> = Traversal<
    NodeData,
    EdgeData,
    IndexType,
    Graph,
    BackwardNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<NodeIndex<IndexType>>,
>;
/// A DFS that treats each directed edge as an undirected edge, i.e. that traverses edge both in forward and backward direction.
pub type UndirectedDfs<NodeData, EdgeData, IndexType, Graph> = Traversal<
    NodeData,
    EdgeData,
    IndexType,
    Graph,
    UndirectedNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<NodeIndex<IndexType>>,
>;

pub struct Traversal<
    NodeData,
    EdgeData,
    IndexType,
    Graph,
    NeighborStrategy,
    QueueStrategy,
    Queue: BidirectedQueue<NodeIndex<IndexType>>,
> {
    queue: Queue,
    order: Vec<IndexType>,
    current_order: IndexType,
    node_data: PhantomData<NodeData>,
    edge_data: PhantomData<EdgeData>,
    graph: PhantomData<Graph>,
    neighbor_strategy: PhantomData<NeighborStrategy>,
    queue_strategy: PhantomData<QueueStrategy>,
}

impl<
        'a,
        NodeData,
        EdgeData,
        IndexType: PrimInt,
        Graph: StaticGraph<IndexType = IndexType>,
        NeighborStrategy: TraversalNeighborStrategy<'a, NodeData, EdgeData, IndexType, Graph>,
        QueueStrategy: TraversalQueueStrategy<IndexType, Queue>,
        Queue: BidirectedQueue<NodeIndex<IndexType>>,
    > Traversal<NodeData, EdgeData, IndexType, Graph, NeighborStrategy, QueueStrategy, Queue>
{
    pub fn new(graph: &Graph, start: NodeIndex<IndexType>) -> Self {
        let mut queue = Queue::default();
        queue.push_back(start);
        let mut order = vec![IndexType::max_value(); graph.node_count()];
        order[start] = IndexType::zero();
        Self {
            queue,
            order,
            current_order: IndexType::one(),
            node_data: Default::default(),
            edge_data: Default::default(),
            graph: Default::default(),
            neighbor_strategy: Default::default(),
            queue_strategy: Default::default(),
        }
    }

    pub fn next(&mut self, graph: &'a Graph) -> Option<NodeIndex<IndexType>> {
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

pub trait TraversalNeighborStrategy<'a, NodeData, EdgeData, IndexType, Graph> {
    type Iterator: Iterator<Item = NodeIndex<IndexType>>;

    fn neighbor_iterator(graph: &'a Graph, node: NodeIndex<IndexType>) -> Self::Iterator;
}

pub trait TraversalQueueStrategy<IndexType, Queue> {
    fn push(queue: &mut Queue, node: NodeIndex<IndexType>);
    fn pop(queue: &mut Queue) -> Option<NodeIndex<IndexType>>;
}

pub struct ForwardNeighborStrategy;
pub type NeighborsIntoNodes<IndexType, Neighbors> = std::iter::Map<
    <Neighbors as IntoIterator>::IntoIter,
    fn(crate::interface::Neighbor<IndexType>) -> NodeIndex<IndexType>,
>;

impl<
        'a,
        NodeData,
        EdgeData,
        IndexType: 'a + PrimInt,
        Graph: NavigableGraph<'a, IndexType = IndexType>,
    > TraversalNeighborStrategy<'a, NodeData, EdgeData, IndexType, Graph>
    for ForwardNeighborStrategy
{
    type Iterator = NeighborsIntoNodes<IndexType, Graph::OutNeighbors>;

    fn neighbor_iterator(graph: &'a Graph, node: NodeIndex<IndexType>) -> Self::Iterator {
        graph
            .out_neighbors(node)
            .unwrap()
            .into_iter()
            .map(|e| e.node_id)
    }
}

pub struct BackwardNeighborStrategy;

impl<
        'a,
        IndexType: 'a + PrimInt,
        NodeData,
        EdgeData,
        Graph: NavigableGraph<'a, IndexType = IndexType>,
    > TraversalNeighborStrategy<'a, NodeData, EdgeData, IndexType, Graph>
    for BackwardNeighborStrategy
{
    type Iterator = NeighborsIntoNodes<IndexType, Graph::InNeighbors>;

    fn neighbor_iterator(graph: &'a Graph, node: NodeIndex<IndexType>) -> Self::Iterator {
        graph
            .in_neighbors(node)
            .unwrap()
            .into_iter()
            .map(|e| e.node_id)
    }
}

pub struct UndirectedNeighborStrategy;
pub type InOutNeighborsIntoNodes<IndexType, OutNeighbors, InNeighbors> = std::iter::Map<
    std::iter::Chain<
        <OutNeighbors as IntoIterator>::IntoIter,
        <InNeighbors as IntoIterator>::IntoIter,
    >,
    fn(crate::interface::Neighbor<IndexType>) -> NodeIndex<IndexType>,
>;

impl<
        'a,
        IndexType: 'a + PrimInt,
        NodeData,
        EdgeData,
        Graph: NavigableGraph<'a, IndexType = IndexType>,
    > TraversalNeighborStrategy<'a, NodeData, EdgeData, IndexType, Graph>
    for UndirectedNeighborStrategy
{
    type Iterator = InOutNeighborsIntoNodes<IndexType, Graph::OutNeighbors, Graph::InNeighbors>;

    fn neighbor_iterator(graph: &'a Graph, node: NodeIndex<IndexType>) -> Self::Iterator {
        graph
            .out_neighbors(node)
            .unwrap()
            .into_iter()
            .chain(graph.in_neighbors(node).unwrap().into_iter())
            .map(|e| e.node_id)
    }
}

pub struct BfsQueueStrategy;

impl<IndexType, Queue: BidirectedQueue<NodeIndex<IndexType>>>
    TraversalQueueStrategy<IndexType, Queue> for BfsQueueStrategy
{
    fn push(queue: &mut Queue, node: NodeIndex<IndexType>) {
        queue.push_back(node)
    }

    fn pop(queue: &mut Queue) -> Option<NodeIndex<IndexType>> {
        queue.pop_front()
    }
}

pub struct DfsQueueStrategy;

impl<IndexType, Queue: BidirectedQueue<NodeIndex<IndexType>>>
    TraversalQueueStrategy<IndexType, Queue> for DfsQueueStrategy
{
    fn push(queue: &mut Queue, node: NodeIndex<IndexType>) {
        queue.push_back(node)
    }

    fn pop(queue: &mut Queue) -> Option<NodeIndex<IndexType>> {
        queue.pop_back()
    }
}
