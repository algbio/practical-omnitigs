//! A crate implementing different omnitig-related algorithms.
#![warn(missing_docs)]

#[macro_use]
extern crate traitsequence_derive;

/// Algorithms to compute the hydrostructure.
pub mod hydrostructure;
/// Algorithms to compute macrotigs.
pub mod macrotigs;
/// Algorithms to extract omnitigs from a graph.
pub mod omnitigs;
/// Algorithms related to restricted reachability queries, like our basic:
/// given an edge e, return everything reachable from the tail of e without using e.
pub mod restricted_reachability;
/// Algorithms to compute unitigs.
pub mod unitigs;

pub use traitgraph;
