/// A subgraph implementation based on bitvectors;
pub mod bit_vector_subgraph;
/// A subgraph implementation that allows to combine multiple subgraphs into one if they are totally ordered by the subset relation.
pub mod incremental_subgraph;
/// A graph implementation based on the `petgraph` crate.
pub mod petgraph_impl;
