//! A crate to represent a bigraph.
//! Bigraphs are graphs where each node is symmetrically mapped to a unique "mirror" node, and each edge is symmetrically mapped to a unique "mirror" edge.
//!
//! Note that bigraphs come in two flavours, node- and edge-centric.
//! A node-centric bigraph has edges that are only distinguished by their endpoints,
//! while an edge-centric bigraph's edges are additionally distinguished by their associated data.
//!
//! This crate implements a simple wrapper around the `traitgraph` crate, adding a vector to represent the node-mirror function.
//! It also implements the edge-mirror function, albeit probably slower than it could be if there was also a vector to map edges.
//!
//! In the context of node-centric genome graphs, nodes usually represent genome strings and a pair of mirrored nodes represent reverse complements of each other.
//! For edge-centric genome graphs, the same holds for the edges.
#![warn(missing_docs)]

/// Abstract algorithms on bigraphs.
pub mod algo;
/// Different implementations of bigraphs.
pub mod implementation;
/// Traits describing the features of bigraphs.
pub mod interface;

pub use traitgraph;
