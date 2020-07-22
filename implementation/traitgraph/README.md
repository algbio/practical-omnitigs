# Traitgraph

[![](http://meritbadge.herokuapp.com/traitgraph)](https://crates.io/crates/traitgraph)
[![](https://docs.rs/traitgraph/badge.svg)](https://docs.rs/traitgraph)
![](https://github.com/algbio/practical-omnitigs/workflows/Tests%20%26%20Lints/badge.svg?branch=master)

A Rust crate to represent and operate on graphs.

The basic principle of this crate is to define all methods on traits, and then implement these for concrete graph representations.
The crate mainly builds on top of [petgraph](https://crates.io/crates/petgraph), and might hopefully be included into petgraph at some point.
