use crate::algo::queue::BidirectedQueue;
use crate::index::{GraphIndex, OptionalGraphIndex};
use crate::interface::NodeOrEdge;
use crate::interface::{GraphBase, ImmutableGraphContainer, NavigableGraph, Neighbor, StaticGraph};
use std::collections::VecDeque;
use std::iter::IntoIterator;
use std::marker::PhantomData;

/// Functions and structures related to univocal traversals.
/// Univocal traversals are traversals along unique out-edges or unique in-edges in a graph.
pub mod univocal_traversal;

/// A normal forward BFS in a directed graph.
pub type PreOrderForwardBfs<'a, Graph> = PreOrderTraversal<
    'a,
    Graph,
    ForwardNeighborStrategy,
    BfsQueueStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;
/// A normal backward BFS in a directed graph.
pub type PreOrderBackwardBfs<'a, Graph> = PreOrderTraversal<
    'a,
    Graph,
    BackwardNeighborStrategy,
    BfsQueueStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;
/// A BFS that treats each directed edge as an undirected edge, i.e. that traverses edge both in forward and backward direction.
pub type PreOrderUndirectedBfs<'a, Graph> = PreOrderTraversal<
    'a,
    Graph,
    UndirectedNeighborStrategy,
    BfsQueueStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;

/// A normal forward DFS in a directed graph.
pub type PreOrderForwardDfs<'a, Graph> = PreOrderTraversal<
    'a,
    Graph,
    ForwardNeighborStrategy,
    DfsQueueStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;
/// A normal backward DFS in a directed graph.
pub type PreOrderBackwardDfs<'a, Graph> = PreOrderTraversal<
    'a,
    Graph,
    BackwardNeighborStrategy,
    DfsQueueStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;
/// A DFS that treats each directed edge as an undirected edge, i.e. that traverses edge both in forward and backward direction.
pub type PreOrderUndirectedDfs<'a, Graph> = PreOrderTraversal<
    'a,
    Graph,
    UndirectedNeighborStrategy,
    DfsQueueStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;

/// A post-order forward DFS in a directed graph.
pub type PostOrderForwardDfs<Graph> = DfsPostOrderTraversal<
    Graph,
    ForwardNeighborStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;
/// A post-order forward DFS in a directed graph.
pub type PostOrderBackwardDfs<Graph> = DfsPostOrderTraversal<
    Graph,
    BackwardNeighborStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;
/// A post-order forward DFS in a directed graph.
pub type PostOrderUndirectedDfs<Graph> = DfsPostOrderTraversal<
    Graph,
    UndirectedNeighborStrategy,
    VecDeque<<Graph as GraphBase>::NodeIndex>,
>;

/// A generic preorder graph traversal.
/// The traversal is generic over the graph implementation,
/// as well as the direction of the search (`NeighborStrategy`),
/// the order of processing (`QueueStrategy`) and the queue implementation itself (`Queue`).
///
/// Moreover, the traversal computes the preorder rank of each visited node.
/// Also, the traversal operates with edge-granularity, meaning that not just nodes are returned by the `next` method, but the traversed edges of each node as well.
/// Additionally, a forbidden subgraph can be passed using the `next_with_forbidden_subgraph` method to disable some edges and nodes in the traversal.
pub struct PreOrderTraversal<
    'a,
    Graph: GraphBase,
    NeighborStrategy: TraversalNeighborStrategy<'a, Graph>,
    QueueStrategy,
    Queue: BidirectedQueue<Graph::NodeIndex>,
> {
    graph: &'a Graph,
    queue: Queue,
    rank: Vec<Graph::OptionalNodeIndex>,
    current_rank: Graph::NodeIndex,
    neighbor_iterator: Option<NeighborStrategy::Iterator>,
    neighbor_strategy: PhantomData<NeighborStrategy>,
    queue_strategy: PhantomData<QueueStrategy>,
}

impl<
        'a,
        Graph: StaticGraph,
        NeighborStrategy: TraversalNeighborStrategy<'a, Graph>,
        QueueStrategy: TraversalQueueStrategy<Graph, Queue>,
        Queue: BidirectedQueue<Graph::NodeIndex>,
    > PreOrderTraversal<'a, Graph, NeighborStrategy, QueueStrategy, Queue>
{
    /// Creates a new traversal that operates on the given graph starting from the given node.
    pub fn new(graph: &'a Graph, start: Graph::NodeIndex) -> Self {
        let mut queue = Queue::default();
        QueueStrategy::push(&mut queue, start);
        let mut rank = vec![Graph::OptionalNodeIndex::new_none(); graph.node_count()];
        rank[start.as_usize()] = Some(0).into();
        Self {
            graph,
            queue,
            rank,
            current_rank: 1.into(),
            neighbor_iterator: None,
            neighbor_strategy: Default::default(),
            queue_strategy: Default::default(),
        }
    }

    /// Creates a new traversal that operates on the given graph.
    /// Does not start the traversal.
    pub fn new_without_start(graph: &'a Graph) -> Self {
        let queue = Queue::default();
        let rank = vec![Graph::OptionalNodeIndex::new_none(); graph.node_count()];
        Self {
            graph,
            queue,
            rank,
            current_rank: 0.into(),
            neighbor_iterator: None,
            neighbor_strategy: Default::default(),
            queue_strategy: Default::default(),
        }
    }

    /// Resets the traversal to start from the given node.
    pub fn reset(&mut self, start: Graph::NodeIndex) {
        self.queue.clear();
        QueueStrategy::push(&mut self.queue, start);
        for rank in &mut self.rank {
            *rank = Graph::OptionalNodeIndex::new_none();
        }
        self.rank[start.as_usize()] = Some(0).into();
        self.current_rank = 1.into();
        self.neighbor_iterator = None;
    }

    /// Resets the traversal to start from the given node without resetting the visited nodes.
    /// Returns the rank of the starting node.
    pub fn continue_traversal_from(&mut self, start: Graph::NodeIndex) -> Graph::NodeIndex {
        debug_assert!(self.queue.is_empty());
        debug_assert!(self.neighbor_iterator.is_none());
        QueueStrategy::push(&mut self.queue, start);
        self.rank[start.as_usize()] = Some(self.current_rank).into();
        let result = self.current_rank;
        self.current_rank = self.current_rank + 1;
        result
    }

    /// Advances the traversal, ignoring all nodes and edges forbidden by `forbidden_subgraph`.
    pub fn next_with_forbidden_subgraph<FN: ForbiddenSubgraph<Graph>>(
        &mut self,
        forbidden_subgraph: &FN,
    ) -> Option<NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>> {
        self.next_internal(forbidden_subgraph)
    }

    #[inline]
    fn next_internal<FS: ForbiddenSubgraph<Graph>>(
        &mut self,
        forbidden_subgraph: &FS,
    ) -> Option<NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>> {
        if let Some(neighbor_iterator) = self.neighbor_iterator.as_mut() {
            for neighbor in neighbor_iterator {
                if forbidden_subgraph.is_edge_forbidden(neighbor.edge_id) {
                    continue;
                }

                if !forbidden_subgraph.is_node_forbidden(neighbor.node_id) {
                    let rank_entry = &mut self.rank[neighbor.node_id.as_usize()];
                    if rank_entry.is_none() {
                        *rank_entry = self.current_rank.into();
                        self.current_rank = self.current_rank + 1;
                        QueueStrategy::push(&mut self.queue, neighbor.node_id);
                    }
                }

                return Some(NodeOrEdge::Edge(neighbor.edge_id));
            }

            self.neighbor_iterator = None;
        }

        if let Some(first) = QueueStrategy::pop(&mut self.queue) {
            debug_assert!(
                !forbidden_subgraph.is_node_forbidden(first),
                "A node became forbidden after being added to the queue. This is not supported."
            );
            self.neighbor_iterator = Some(NeighborStrategy::neighbor_iterator(self.graph, first));

            Some(NodeOrEdge::Node(first))
        } else {
            None
        }
    }

    /// Returns the rank of the given node, or `None` if the node has not yet been visited.
    pub fn rank_of(&self, node: Graph::NodeIndex) -> Option<Graph::NodeIndex> {
        let rank = self.rank[node.as_usize()];
        rank.into()
    }
}
impl<
        'a,
        Graph: StaticGraph,
        NeighborStrategy: TraversalNeighborStrategy<'a, Graph>,
        QueueStrategy: TraversalQueueStrategy<Graph, Queue>,
        Queue: BidirectedQueue<Graph::NodeIndex>,
    > Iterator for PreOrderTraversal<'a, Graph, NeighborStrategy, QueueStrategy, Queue>
{
    type Item = NodeOrEdge<Graph::NodeIndex, Graph::EdgeIndex>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_internal(&NoForbiddenSubgraph)
    }
}

/// A generic depth first postorder graph traversal.
/// The traversal is generic over the graph implementation,
/// as well as the direction of the search (`NeighborStrategy`)
/// and the queue implementation (`Queue`).
///
/// Moreover, the traversal computes the postorder rank of each visited node.
/// This traversal operates with node-granularity, meaning that the `next` method returns nodes.
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
    /// Creates a new traversal that operates on the given graph, starting from the given node.
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

    /// Creates a new traversal that operates on the given graph.
    /// There is no starting node given, and to start the search, one of the `reset` methods needs to be used.
    pub fn new_without_start(graph: &Graph) -> Self {
        let queue = Queue::default();
        let rank = vec![Graph::OptionalNodeIndex::new_none(); graph.node_count()];
        Self {
            queue,
            rank,
            current_rank: 0.into(),
            graph: Default::default(),
            neighbor_strategy: Default::default(),
        }
    }

    /// Resets the traversal to start from the given node.
    pub fn reset(&mut self, start: Graph::NodeIndex) {
        self.queue.clear();
        self.queue.push_back(start);
        for rank in &mut self.rank {
            *rank = Graph::OptionalNodeIndex::new_none();
        }
        self.current_rank = 0.into();
    }

    /// Resets the traversal to start from the given node without resetting the visited nodes.
    pub fn continue_traversal_from(&mut self, start: Graph::NodeIndex) {
        assert!(self.queue.is_empty());
        self.queue.push_back(start);
    }

    /// Computes and returns the next node in depth-first search postorder.
    pub fn next(&mut self, graph: &'a Graph) -> Option<Graph::NodeIndex> {
        while let Some(first) = self.queue.pop_back() {
            let rank_entry = &mut self.rank[first.as_usize()];
            if *rank_entry == Self::explored_rank() {
                debug_assert_ne!(self.current_rank.into(), Self::explored_rank());
                *rank_entry = self.current_rank.into();
                self.current_rank = self.current_rank + 1;

                return Some(first);
            } else if rank_entry.is_none() {
                self.queue.push_back(first);
                *rank_entry = Self::explored_rank();

                for neighbor in NeighborStrategy::neighbor_iterator(graph, first) {
                    let rank_entry = &mut self.rank[neighbor.node_id.as_usize()];
                    if rank_entry.is_none() {
                        self.queue.push_back(neighbor.node_id);
                    }
                }
            }
        }

        None
    }

    /// Returns the rank of a node in depth-first search postorder, or `None` if the node has not yet been processed completely.
    pub fn rank_of(&self, node: Graph::NodeIndex) -> Option<Graph::NodeIndex> {
        let rank = self.rank[node.as_usize()];
        rank.into()
    }

    fn explored_rank() -> Graph::OptionalNodeIndex {
        Some(Graph::OptionalNodeIndex::new_none().as_usize_unchecked() - 1).into()
    }
}

/// A type with this trait can tell if a node or edge is forbidden in a graph traversal.
pub trait ForbiddenSubgraph<Graph: GraphBase> {
    /// Returns true if the given node is forbidden.
    fn is_node_forbidden(&self, node: Graph::NodeIndex) -> bool;

    /// Returns true if the given edge is forbidden.
    fn is_edge_forbidden(&self, edge: Graph::EdgeIndex) -> bool;
}

/// A type that defines the strategy for computing the neighborhood of a node or edge, i.e. forward, backward or undirected.
pub trait TraversalNeighborStrategy<'a, Graph: GraphBase> {
    /// The iterator type used to iterate over the neighbors of a node.
    type Iterator: Iterator<Item = Neighbor<Graph::NodeIndex, Graph::EdgeIndex>>;
    /// The iterator type used to iterate over the neighbors of an edge.
    type EdgeNeighborIterator: Iterator<Item = Graph::NodeIndex>;

    /// Returns an iterator over the neighbors of a given node.
    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator;

    /// Returns an iterator over the neighbors of an edge.
    fn edge_neighbor_iterator(
        graph: &'a Graph,
        edge: Graph::EdgeIndex,
    ) -> Self::EdgeNeighborIterator;
}

/// A type that defines the order of node processing in a traversal, i.e. queue-based or stack-based.
pub trait TraversalQueueStrategy<Graph: GraphBase, Queue: BidirectedQueue<Graph::NodeIndex>> {
    /// Insert a node into the queue.
    fn push(queue: &mut Queue, node: Graph::NodeIndex);
    /// Remove and return a node from the queue.
    fn pop(queue: &mut Queue) -> Option<Graph::NodeIndex>;
}

/// A type implementing [ForbiddenSubgraph](ForbiddenSubgraph) that allows all nodes in a graph traversal.
pub struct NoForbiddenSubgraph;
impl<Graph: GraphBase> ForbiddenSubgraph<Graph> for NoForbiddenSubgraph {
    fn is_node_forbidden(&self, _: Graph::NodeIndex) -> bool {
        false
    }

    fn is_edge_forbidden(&self, _: Graph::EdgeIndex) -> bool {
        false
    }
}

/// A type implementing [ForbiddenSubgraph](ForbiddenSubgraph) that allows all nodes set to true in a boolean vector.
pub struct AllowedNodesForbiddenSubgraph<'a> {
    allowed_nodes: &'a [bool],
}
impl<'a> AllowedNodesForbiddenSubgraph<'a> {
    /// Creates a new `AllowedNodesForbiddenSubgraph` with the given boolean vector that contains `true` for each allowed node and `false` for each forbidden node.
    pub fn new(allowed_nodes: &'a [bool]) -> Self {
        Self { allowed_nodes }
    }
}
impl<'a, Graph: GraphBase> ForbiddenSubgraph<Graph> for AllowedNodesForbiddenSubgraph<'a> {
    fn is_node_forbidden(&self, node: Graph::NodeIndex) -> bool {
        !self.allowed_nodes[node.as_usize()]
    }

    fn is_edge_forbidden(&self, _: Graph::EdgeIndex) -> bool {
        false
    }
}

/// A [ForbiddenSubgraph](ForbiddenSubgraph) that forbids a single edge.
pub struct ForbiddenEdge<EdgeIndex> {
    edge_id: EdgeIndex,
}
impl<EdgeIndex> ForbiddenEdge<EdgeIndex> {
    /// Construct a new `ForbiddenEdge` that forbids the given edge.
    pub fn new(forbidden_edge: EdgeIndex) -> Self {
        Self {
            edge_id: forbidden_edge,
        }
    }
}
impl<Graph: GraphBase> ForbiddenSubgraph<Graph> for ForbiddenEdge<Graph::EdgeIndex> {
    fn is_node_forbidden(&self, _: Graph::NodeIndex) -> bool {
        false
    }

    fn is_edge_forbidden(&self, edge: Graph::EdgeIndex) -> bool {
        edge == self.edge_id
    }
}

/// A [ForbiddenSubgraph](ForbiddenSubgraph) that forbids a single node.
pub struct ForbiddenNode<NodeIndex> {
    node_id: NodeIndex,
}
impl<NodeIndex> ForbiddenNode<NodeIndex> {
    /// Construct a new `ForbiddenNode` that forbids the given node.
    pub fn new(forbidden_node: NodeIndex) -> Self {
        Self {
            node_id: forbidden_node,
        }
    }
}
impl<Graph: GraphBase> ForbiddenSubgraph<Graph> for ForbiddenNode<Graph::NodeIndex> {
    fn is_node_forbidden(&self, node: Graph::NodeIndex) -> bool {
        node == self.node_id
    }

    fn is_edge_forbidden(&self, _: Graph::EdgeIndex) -> bool {
        false
    }
}

/// A neighbor strategy that traverses all outgoing edges of a node.
pub struct ForwardNeighborStrategy;
/*pub type NeighborsIntoNodes<NodeIndex, EdgeIndex, Neighbors> = std::iter::Map<
    <Neighbors as IntoIterator>::IntoIter,
    fn(crate::interface::Neighbor<NodeIndex, EdgeIndex>) -> NodeIndex,
>;*/

impl<'a, Graph: NavigableGraph<'a> + ImmutableGraphContainer> TraversalNeighborStrategy<'a, Graph>
    for ForwardNeighborStrategy
{
    type Iterator = Graph::OutNeighbors;
    type EdgeNeighborIterator = std::iter::Once<Graph::NodeIndex>;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator {
        graph.out_neighbors(node)
    }

    fn edge_neighbor_iterator(
        graph: &'a Graph,
        edge: Graph::EdgeIndex,
    ) -> Self::EdgeNeighborIterator {
        std::iter::once(graph.edge_endpoints(edge).to_node)
    }
}

/// A neighbor strategy that traverses all incoming edges of a node.
pub struct BackwardNeighborStrategy;

impl<'a, Graph: NavigableGraph<'a> + ImmutableGraphContainer> TraversalNeighborStrategy<'a, Graph>
    for BackwardNeighborStrategy
{
    type Iterator = Graph::InNeighbors;
    type EdgeNeighborIterator = std::iter::Once<Graph::NodeIndex>;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator {
        graph.in_neighbors(node)
    }

    fn edge_neighbor_iterator(
        graph: &'a Graph,
        edge: Graph::EdgeIndex,
    ) -> Self::EdgeNeighborIterator {
        std::iter::once(graph.edge_endpoints(edge).from_node)
    }
}

/// A neighbor strategy that traverses all incoming and all outgoing edges of a node.
pub struct UndirectedNeighborStrategy;
type InOutNeighborsChain<OutNeighbors, InNeighbors> = std::iter::Chain<
    <OutNeighbors as IntoIterator>::IntoIter,
    <InNeighbors as IntoIterator>::IntoIter,
>;

impl<'a, Graph: NavigableGraph<'a> + ImmutableGraphContainer> TraversalNeighborStrategy<'a, Graph>
    for UndirectedNeighborStrategy
{
    type Iterator = InOutNeighborsChain<Graph::OutNeighbors, Graph::InNeighbors>;
    type EdgeNeighborIterator =
        std::iter::Chain<std::iter::Once<Graph::NodeIndex>, std::iter::Once<Graph::NodeIndex>>;

    fn neighbor_iterator(graph: &'a Graph, node: Graph::NodeIndex) -> Self::Iterator {
        graph.out_neighbors(node).chain(graph.in_neighbors(node))
    }

    fn edge_neighbor_iterator(
        graph: &'a Graph,
        edge: Graph::EdgeIndex,
    ) -> Self::EdgeNeighborIterator {
        std::iter::once(graph.edge_endpoints(edge).to_node)
            .chain(std::iter::once(graph.edge_endpoints(edge).from_node))
    }
}

/// A queue strategy that works by the first-in first-out principle.
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

/// A queue strategy that works by the last-in first-out principle.
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
    use std::collections::VecDeque;

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
            DfsPostOrderTraversal::<_, ForwardNeighborStrategy, VecDeque<_>>::new(&graph, n0);
        debug_assert_eq!(
            graph.out_neighbors(n0).map(|n| n.node_id).next(),
            Some(3.into())
        );
        debug_assert_eq!(ordering.next(&graph), Some(n3));
        debug_assert_eq!(ordering.next(&graph), Some(n2));
        debug_assert_eq!(ordering.next(&graph), Some(n1));
        debug_assert_eq!(ordering.next(&graph), Some(n0));
        debug_assert_eq!(ordering.next(&graph), None);
    }
}
