#![warn(missing_docs)]
//! This crate offers traits for abstract graph algorithms as well as implementations of these traits.

pub use traitsequence;

/// Abstract graph algorithms.
pub mod algo;
/// Different implementations of the graph traits.
pub mod implementation;
/// Traits and a default implementation for graph indices.
pub mod index;
pub mod interface;
/// Methods for reading and writing graphs.
pub mod io;
/// Traits and implementations of node- and edge-centric walks.
pub mod walks;
