//! A crate implementing different omnitig-related algorithms.
#![warn(missing_docs)]

/// Algorithms to compute macrotigs.
pub mod macrotigs;
/// Algorithms related to restricted reachability queries, like our basic:
/// given an edge e, return everything reachable from the tail of e without using e.
pub mod restricted_reachability;
/// Algorithms to compute unitigs.
pub mod unitigs;

pub use traitgraph;
