//! A crate implementing different omnitig-related algorithms.
#![warn(missing_docs)]

#[macro_use]
extern crate log;

/// Preprocess a hamiltonian circuit problem using node-centric omnitigs.
pub mod hamiltonian;
/// Algorithms to compute the hydrostructure.
pub mod hydrostructure;
/// Algorithms to compute macrotigs.
pub mod macrotigs;
/// Algorithms to compute maximal safe walks under the node-covering node-visible 1-circular walk model.
pub mod node_covering_node_visible_one_circular_safe;
/// Algorithms to extract omnitigs from a graph.
pub mod omnitigs;
/// Algorithms related to restricted reachability queries, like our basic:
/// given an edge e, return everything reachable from the tail of e without using e.
pub mod restricted_reachability;
/// Algorithms to compute unitigs.
pub mod unitigs;
/// Core functions on walks for omnitig-like algorithms.
pub mod walks;

pub use traitgraph;
