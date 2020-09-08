//! This crate abstracts over the representation of a genome string, allowing for different implementations that are catered to different use-cases.
#![warn(missing_docs)]

/// Different implementations of genome string representations.
pub mod implementation;
pub mod interface;
