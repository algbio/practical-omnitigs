use crate::algo::queue::BidirectedQueue;
use crate::index::{GraphIndex, OptionalGraphIndex};
use crate::interface::{GraphBase, NavigableGraph, StaticGraph};
use std::collections::LinkedList;
use std::iter::IntoIterator;
use std::marker::PhantomData;

/// A normal forward BFS in a directed graph.
pub type PreOrderForwardBfs<Graph> = PreOrderTraversal<
    Graph,
    ForwardNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A normal backward BFS in a directed graph.
pub type PreOrderBackwardBfs<Graph> = PreOrderTraversal<
    Graph,
    BackwardNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A BFS that treats each directed edge as an undirected edge, i.e. that traverses edge both in forward and backward direction.
pub type PreOrderUndirectedBfs<Graph> = PreOrderTraversal<
    Graph,
    UndirectedNeighborStrategy,
    BfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;

/// A normal forward DFS in a directed graph.
pub type PreOrderForwardDfs<Graph> = PreOrderTraversal<
    Graph,
    ForwardNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A normal backward DFS in a directed graph.
pub type PreOrderBackwardDfs<Graph> = PreOrderTraversal<
    Graph,
    BackwardNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A DFS that treats each directed edge as an undirected edge, i.e. that traverses edge both in forward and backward direction.
pub type PreOrderUndirectedDfs<Graph> = PreOrderTraversal<
    Graph,
    UndirectedNeighborStrategy,
    DfsQueueStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;

/// A post-order forward DFS in a directed graph.
pub type PostOrderForwardDfs<Graph> = DfsPostOrderTraversal<
    Graph,
    ForwardNeighborStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A post-order forward DFS in a directed graph.
pub type PostOrderBackwardDfs<Graph> = DfsPostOrderTraversal<
    Graph,
    BackwardNeighborStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;
/// A post-order forward DFS in a directed graph.
pub type PostOrderUndirectedDfs<Graph> = DfsPostOrderTraversal<
    Graph,
    UndirectedNeighborStrategy,
    LinkedList<<Graph as GraphBase>::NodeIndex>,
>;

pub struct PreOrderTraversal<
    Graph: GraphBase,
    NeighborStrategy,
    QueueStrategy,
    Queue: BidirectedQueue<Graph::NodeIndex>,
> {
    queue: Queue,
    rank: Vec<Graph::OptionalNodeIndex>,
    current_rank: Graph::NodeIndex,
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
    > PreOrderTraversal<Graph, NeighborStrategy, QueueStrategy, Queue>
{
    pub fn new(graph: &Graph, start: Graph::NodeIndex) -> Self {
        let mut queue = Queue::default();
        QueueStrategy::push(&mut queue, start);
        let mut rank = vec![Graph::OptionalNodeIndex::new_none(); graph.node_count()];
        rank[start.as_usize()] = Some(0).into();
        Self {
            queue,
            rank,
            current_rank: 1.into(),
            graph: Default::default(),
            neighbor_strategy: Default::default(),
            queue_strategy: Default::default(),
        }
    }

    pub fn next(&mut self, graph: &'a Graph) -> Option<Graph::NodeIndex> {
        self.next_internal(graph, &NoForbiddenNodes)
    }

    pub fn next_with_forbidden_nodes<FN: ForbiddenNodes<Graph>>(
        &mut self,
        graph: &'a Graph,
        forbidden_nodes: &FN,
    ) -> Option<Graph::NodeIndex> {
        self.next_internal(graph, forbidden_nodes)
    }

    #[inline]
    fn next_internal<FN: ForbiddenNodes<Graph>>(
        &mut self,
        graph: &'a Graph,
        forbidden_nodes: &FN,
    ) -> Option<Graph::NodeIndex> {
        if let Some(first) = QueueStrategy::pop(&mut self.queue) {
            assert!(
                !forbidden_nodes.is_forbidden(first),
                "A node became forbidden after being added to the queue. This is not supported."
            );

            for neighbor in NeighborStrategy::neighbor_iterator(graph, first) {
                if forbidden_nodes.is_forbidden(neighbor) {
                    continue;
                }

                let rank_entry = &mut self.rank[neighbor.as_usize()];
                if *rank_entry == None.into() {
                    *rank_entry = self.current_rank.into();
                    self.current_rank = self.current_rank + 1;
                    QueueStrategy::push(&mut self.queue, neighbor);
                }
            }

            Some(first)
        } else {
            None
        }
    }

    pub fn rank_of(&self, node: Graph::NodeIndex) -> Option<Graph::NodeIndex> {
        let rank = self.rank[node.as_usize()];
        rank.into()
    }
}

pub struct DfsPostOrderTraversal<
    Graph: GraphBase,
    NeighborStrategy,
    Queue: BidirectedQueue<Graph::NodeIndex>,
> {
    queue: Queue,
    rank: Vec<Graph::OptionalNodeIndex>,
    current_rank: Graph::NodeIndex,
    graph: PhantomData<Graph>,
    neighbor_strategy: PhantomData<NeighborStrategy>,
}

impl<
        'a,
        Graph: StaticGraph,
        NeighborStrategy: TraversalNeighborStrategy<'a, Graph>,
        Queue: BidirectedQueue<Graph::NodeIndex>,
    > DfsPostOrderTraversal<Graph, NeighborStrategy, Queue>
{
    pub fn new(graph: &Graph, start: Graph::NodeIndex) -> Self {
        let mut queue = Queue::default();
        queue.push_back(start);
        let rank = vec![Graph::OptionalNodeIndex::new_none(); graph.node_count()];
        Self {
            queue,
            rank,
            current_rank: 0.into(),
            graph: Default::default(),
            neighbor_strategy: Default::default(),
        }
    }

    pub fn reset(&mut self, start: Graph::NodeIndex) {
        self.queue.clear();
        self.queue.push_back(start);
        for rank in &mut self.rank {
            *rank = Graph::OptionalNodeIndex::new_none();
        }
        self.current_rank = 0.into();
    }

    pub fn next(&mut self, graph: &'a Graph) -> Option<Graph::NodeIndex> {
        while let Some(first) = self.queue.pop_back() {
            let rank_entry = &mut self.rank[first.as_usize()];
            if *rank_entry == Self::explored_rank() {
                assert_ne!(self.current_rank.into(), Self::explored_rank());
                *rank_entry = self.current_rank.into();
                self.current_rank = self.current_rank + 1;

                return Some(first);
            } else if *rank_entry == None.into() {
                self.queue.push_back(first);
                *rank_entry = Self::explored_rank();

                for neighbor in NeighborStrategy::neighbor_iterator(graph, first) {
                    let rank_entry = &mut self.rank[neighbor.as_usize()];
                    if *rank_entry == None.into() {
                        self.queue.push_back(neighbor);
                    }
                }
            }
        }

        None
    }

    pub fn rank_of(&self, node: Graph::NodeIndex) -> Option<Graph::NodeIndex> {
        let rank = self.rank[node.as_usize()];
        rank.into()
    }

    fn explored_rank() -> Graph::OptionalNodeIndex {
        Some(Graph::OptionalNodeIndex::new_none().as_usize_unchecked() - 1).into()
    }
}

/// A type with this trait can tell if a node is forbidden in a graph traversal.
pub trait ForbiddenNodes<Graph: GraphBase> {
    fn is_forbidden(&self, node: Graph::NodeIndex) -> bool;
}

pub trait TraversalNeighborStrategy<'a, Graph: GraphBase> {
    type Iterator: Iterator<Item = <Graph as GraphBase>::NodeIndex>;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator;
}

pub trait TraversalQueueStrategy<Graph: GraphBase, Queue: BidirectedQueue<Graph::NodeIndex>> {
    fn push(queue: &mut Queue, node: Graph::NodeIndex);
    fn pop(queue: &mut Queue) -> Option<Graph::NodeIndex>;
}

/// A type implementing [ForbiddenNodes](ForbiddenNodes) that allows all nodes in a graph traversal.
pub struct NoForbiddenNodes;
impl<Graph: GraphBase> ForbiddenNodes<Graph> for NoForbiddenNodes {
    fn is_forbidden(&self, _: Graph::NodeIndex) -> bool {
        false
    }
}

/// A type implementing [ForbiddenNodes](ForbiddenNodes) that allows all nodes set to true in a boolean vector.
pub struct AllowedForbiddenNodes<'a> {
    allowed_nodes: &'a [bool],
}
impl<'a> AllowedForbiddenNodes<'a> {
    pub fn new(allowed_nodes: &'a [bool]) -> Self {
        Self { allowed_nodes }
    }
}
impl<'a, Graph: GraphBase> ForbiddenNodes<Graph> for AllowedForbiddenNodes<'a> {
    fn is_forbidden(&self, node: Graph::NodeIndex) -> bool {
        !self.allowed_nodes[node.as_usize()]
    }
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
        graph.out_neighbors(node).map(|e| e.node_id)
    }
}

pub struct BackwardNeighborStrategy;

impl<'a, Graph: NavigableGraph<'a>> TraversalNeighborStrategy<'a, Graph>
    for BackwardNeighborStrategy
{
    type Iterator = NeighborsIntoNodes<Graph::NodeIndex, Graph::EdgeIndex, Graph::InNeighbors>;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator {
        graph.in_neighbors(node).map(|e| e.node_id)
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
            .chain(graph.in_neighbors(node))
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

#[cfg(test)]
mod test {
    use crate::algo::traversal::{DfsPostOrderTraversal, ForwardNeighborStrategy};
    use crate::implementation::petgraph_impl;
    use crate::interface::{MutableGraphContainer, NavigableGraph};
    use std::collections::LinkedList;

    #[test]
    fn test_postorder_traversal_simple() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(0);
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        graph.add_edge(n0, n1, 10);
        graph.add_edge(n1, n2, 11);
        graph.add_edge(n2, n3, 12);
        graph.add_edge(n3, n0, 13);
        graph.add_edge(n1, n0, 14);
        graph.add_edge(n2, n1, 15);
        graph.add_edge(n3, n2, 16);
        graph.add_edge(n0, n3, 17);

        let mut ordering =
            DfsPostOrderTraversal::<_, ForwardNeighborStrategy, LinkedList<_>>::new(&graph, n0);
        assert_eq!(
            graph.out_neighbors(n0).map(|n| n.node_id).next(),
            Some(3.into())
        );
        assert_eq!(ordering.next(&graph), Some(n3));
        assert_eq!(ordering.next(&graph), Some(n2));
        assert_eq!(ordering.next(&graph), Some(n1));
        assert_eq!(ordering.next(&graph), Some(n0));
        assert_eq!(ordering.next(&graph), None);
    }
}
