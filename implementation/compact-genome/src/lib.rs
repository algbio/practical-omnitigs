//! This crate abstracts over the representation of a genome string, allowing for different implementations that are catered to different use-cases.
#![warn(missing_docs)]

/// Different implementations of genome string representations.
pub mod implementation;
pub mod interface;

const ASCII_A: u8 = b'A';
const ASCII_C: u8 = b'C';
const ASCII_G: u8 = b'G';
const ASCII_T: u8 = b'T';
