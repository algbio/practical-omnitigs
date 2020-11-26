//! A crate to represent genome graphs.
//!
//! Genome graphs are typically node- or edge-centric bigraphs that store genome strings on their nodes or edges respectively.
//! This crate offers type aliases using the `bigraph` crate to easily define genome graphs, as well as methods for reading and writing them.
//!
//! Currently, the format for input and output is the [bcalm2 fasta format](https://github.com/GATB/bcalm).
#![warn(missing_docs)]
#![recursion_limit = "1024"]
#[macro_use]
extern crate error_chain;
#[macro_use]
extern crate log;

/// Contains the error types used by this crate.
pub mod error;
/// Contains functions for reading and writing genome graphs.
pub mod io;
/// Contains type aliases for genome graphs.
pub mod types;

pub use bigraph;
pub use compact_genome;
