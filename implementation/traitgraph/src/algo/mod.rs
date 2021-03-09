/// Algorithms related to graph components, i.e. finding the strongly or weakly connected components of a graph or checking if a graph is strongly connected.
pub mod components;
/// Dijkstra's shortest path algorithm.
pub mod dijkstra;
/// Algorithms related to Eulerian graphs.
pub mod eulerian;
/// Algorithms to create certain parameterisable graph classes, like binary trees.
pub mod predefined_graphs;
/// A trait for bidirected queues to abstract over the different implementations in the standard library.
pub mod queue;
/// Algorithms for graph traversals, i.e. preorder breadth or depth first search as well as postorder depth first search.
pub mod traversal;
