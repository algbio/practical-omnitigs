use crate::algo::queue::BidirectedQueue;
use crate::{GraphBase, GraphIndex, NavigableGraph, OptionalGraphIndex, StaticGraph};
use std::collections::LinkedList;
use std::marker::PhantomData;

/// A normal forward BFS in a directed graph.
pub type ForwardBfs<Graph> = Traversal<
    Graph,
    ForwardNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A normal backward BFS in a directed graph.
pub type BackwardBfs<Graph> = Traversal<
    Graph,
    BackwardNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A BFS that treats each directed edge as an undirected edge, i.e. that traverses edge both in forward and backward direction.
pub type UndirectedBfs<Graph> = Traversal<
    Graph,
    UndirectedNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;

/// A normal forward DFS in a directed graph.
pub type ForwardDfs<Graph> = Traversal<
    Graph,
    ForwardNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A normal backward DFS in a directed graph.
pub type BackwardDfs<Graph> = Traversal<
    Graph,
    BackwardNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A DFS that treats each directed edge as an undirected edge, i.e. that traverses edge both in forward and backward direction.
pub type UndirectedDfs<Graph> = Traversal<
    Graph,
    UndirectedNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;

pub struct Traversal<
    Graph: GraphBase,
    NeighborStrategy,
    QueueStrategy,
    Queue: BidirectedQueue<Graph::NodeIndex>,
> {
    queue: Queue,
    order: Vec<Graph::OptionalNodeIndex>,
    current_order: Graph::NodeIndex,
    graph: PhantomData<Graph>,
    neighbor_strategy: PhantomData<NeighborStrategy>,
    queue_strategy: PhantomData<QueueStrategy>,
}

impl<
        'a,
        Graph: StaticGraph,
        NeighborStrategy: TraversalNeighborStrategy<'a, Graph>,
        QueueStrategy: TraversalQueueStrategy<Graph, Queue>,
        Queue: BidirectedQueue<Graph::NodeIndex>,
    > Traversal<Graph, NeighborStrategy, QueueStrategy, Queue>
{
    pub fn new(graph: &Graph, start: Graph::NodeIndex) -> Self {
        let mut queue = Queue::default();
        queue.push_back(start);
        let mut order = vec![Graph::OptionalNodeIndex::new_none(); graph.node_count()];
        order[start.as_usize()] = Some(0).into();
        Self {
            queue,
            order,
            current_order: 1.into(),
            graph: Default::default(),
            neighbor_strategy: Default::default(),
            queue_strategy: Default::default(),
        }
    }

    pub fn next(&mut self, graph: &'a Graph) -> Option<Graph::NodeIndex> {
        if let Some(first) = self.queue.pop_front() {
            for neighbor in graph.out_neighbors(first) {
                let order_entry = &mut self.order[neighbor.node_id.as_usize()];
                if *order_entry == None.into() {
                    *order_entry = self.current_order.into();
                    self.current_order = self.current_order + 1;
                    self.queue.push_back(neighbor.node_id);
                }
            }

            Some(first)
        } else {
            None
        }
    }

    pub fn order_of(&self, node: Graph::NodeIndex) -> Option<Graph::NodeIndex> {
        let order = self.order[node.as_usize()];
        order.into()
    }
}

pub trait TraversalNeighborStrategy<'a, Graph: GraphBase> {
    type Iterator: Iterator<Item = <Graph as GraphBase>::NodeIndex>;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator;
}

pub trait TraversalQueueStrategy<Graph: GraphBase, Queue: BidirectedQueue<Graph::NodeIndex>> {
    fn push(queue: &mut Queue, node: Graph::NodeIndex);
    fn pop(queue: &mut Queue) -> Option<Graph::NodeIndex>;
}

pub struct ForwardNeighborStrategy;
pub type NeighborsIntoNodes<NodeIndex, EdgeIndex, Neighbors> = std::iter::Map<
    <Neighbors as IntoIterator>::IntoIter,
    fn(crate::interface::Neighbor<NodeIndex, EdgeIndex>) -> NodeIndex,
>;

impl<'a, Graph: NavigableGraph<'a>> TraversalNeighborStrategy<'a, Graph>
    for ForwardNeighborStrategy
{
    type Iterator = NeighborsIntoNodes<Graph::NodeIndex, Graph::EdgeIndex, Graph::OutNeighbors>;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator {
        graph.out_neighbors(node).into_iter().map(|e| e.node_id)
    }
}

pub struct BackwardNeighborStrategy;

impl<'a, Graph: NavigableGraph<'a>> TraversalNeighborStrategy<'a, Graph>
    for BackwardNeighborStrategy
{
    type Iterator = NeighborsIntoNodes<Graph::NodeIndex, Graph::EdgeIndex, Graph::InNeighbors>;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator {
        graph.in_neighbors(node).into_iter().map(|e| e.node_id)
    }
}

pub struct UndirectedNeighborStrategy;
pub type InOutNeighborsIntoNodes<NodeIndex, EdgeIndex, OutNeighbors, InNeighbors> = std::iter::Map<
    std::iter::Chain<
        <OutNeighbors as IntoIterator>::IntoIter,
        <InNeighbors as IntoIterator>::IntoIter,
    >,
    fn(crate::interface::Neighbor<NodeIndex, EdgeIndex>) -> NodeIndex,
>;

impl<'a, Graph: NavigableGraph<'a>> TraversalNeighborStrategy<'a, Graph>
    for UndirectedNeighborStrategy
{
    type Iterator = InOutNeighborsIntoNodes<
        Graph::NodeIndex,
        Graph::EdgeIndex,
        Graph::OutNeighbors,
        Graph::InNeighbors,
    >;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator {
        graph
            .out_neighbors(node)
            .into_iter()
            .chain(graph.in_neighbors(node).into_iter())
            .map(|e| e.node_id)
    }
}

pub struct BfsQueueStrategy;

impl<Graph: GraphBase, Queue: BidirectedQueue<Graph::NodeIndex>>
    TraversalQueueStrategy<Graph, Queue> for BfsQueueStrategy
{
    fn push(queue: &mut Queue, node: Graph::NodeIndex) {
        queue.push_back(node)
    }

    fn pop(queue: &mut Queue) -> Option<Graph::NodeIndex> {
        queue.pop_front()
    }
}

pub struct DfsQueueStrategy;

impl<Graph: GraphBase, Queue: BidirectedQueue<Graph::NodeIndex>>
    TraversalQueueStrategy<Graph, Queue> for DfsQueueStrategy
{
    fn push(queue: &mut Queue, node: Graph::NodeIndex) {
        queue.push_back(node)
    }

    fn pop(queue: &mut Queue) -> Option<Graph::NodeIndex> {
        queue.pop_back()
    }
}
