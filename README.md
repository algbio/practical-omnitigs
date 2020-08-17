# Practical [Omnitigs](https://www.liebertpub.com/doi/full/10.1089/cmb.2016.0141)

![](https://github.com/algbio/practical-omnitigs/workflows/Tests%20%26%20Lints/badge.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/algbio/practical-omnitigs/badge.svg?branch=master)](https://coveralls.io/github/algbio/practical-omnitigs?branch=master)

This is a repository to conduct experiments with omnitig-related models in the context of genome assembly.
The algorithms are implemented in *Rust*, and around that we built a *snakemake* toolchain to conduct experiments.
To ensure reproducibility, we wrapped everything into a *conda* environment.

## Usage

### Required Software

 * `conda >= 4.8.3` (lower might be possible, but has not been tested)

### Setup

First, set up the conda environment of this project.
```commandline
conda env create -f environment.yml
```

Then, activate the environment.
```commandline
source activate practical-omnitigs
```

### Running Experiments

Make sure that you are in the right conda environment (should be `practical-omnitigs`).
```commandline
conda info
```

Subsequently, experiments can be run using `snakemake`.
```commandline
snakemake --cores all <experiment>
```

Valid experiments are:
 * `selftest`: Check if you have conda set up correctly. It prints the version of `snakemake` `conda` and `wget`. The versions of `snakemake` and `conda` should match the definition in `/environment.yaml` and the version of `wget` should match the definition in `/config/conda-selftest-env.yaml`.
 * `test`: execute all integration tests of this project on a single small sample genome.
 * `test_all`: execute all integration tests of this project on all defined genomes (potentially very large).

The experiments are run inside a conda environment that is set up by snakemake.
This ensures reproducibility of the results and automates the installation of required tools.

### Using the Implementation Directly

The Rust code written for this project includes a command line interface that can be used directly.
For documentation on how to use it, please refer to the documentation of the [cli crate][cli crate].

## Technical Information

### Directory Structure

 * `.github`: GitHub workflows for continuous testing.
 * `.idea`: Configuration for the IntelliJ IDEA integrated development environment.
 * `config`: All config files related to the experiments, including conda environments and experiment declarations.
 * `data`: Data used and produced by the experiments.
 * `implementation`: The algorithms that we are testing. Everything is written in Rust, so the folder contains a Rust workspace. 

### Implementation

The algorithms of this project are implemented in Rust.
We split the implementation into multiple library crates to increase the reusability of our code.
On top of that, the `cli` crate provides all implemented functionality via a command line interface.
Refer to [its documentation][cli crate] for more information.

This project contains the following crates.

 * [bigraph](https://github.com/algbio/practical-omnitigs/tree/master/implementation/bigraph)
   [![](http://meritbadge.herokuapp.com/bigraph)](https://crates.io/crates/bigraph)
   [![](https://docs.rs/bigraph/badge.svg)](https://docs.rs/bigraph)
   Representation of and operations on bigraphs.
 * [cli][cli crate]
   The command line interface of our implementation.
 * [compact-genome](https://github.com/algbio/practical-omnitigs/tree/master/implementation/compact-genome)
   [![](http://meritbadge.herokuapp.com/compact-genome)](https://crates.io/crates/compact-genome)
   [![](https://docs.rs/compact-genome/badge.svg)](https://docs.rs/compact-genome)
   Representation of and operations on genome strings.
 * [genome-graph](https://github.com/algbio/practical-omnitigs/tree/master/implementation/genome-graph)
   [![](http://meritbadge.herokuapp.com/genome-graph)](https://crates.io/crates/genome-graph)
   [![](https://docs.rs/genome-graph/badge.svg)](https://docs.rs/genome-graph)
   Representation of and operations on genome graphs.
 * [traitgraph](https://github.com/algbio/practical-omnitigs/tree/master/implementation/traitgraph)
   [![](http://meritbadge.herokuapp.com/traitgraph)](https://crates.io/crates/traitgraph)
   [![](https://docs.rs/traitgraph/badge.svg)](https://docs.rs/traitgraph)
   Abstract algorithms on graphs, independent of the graph representation.

Except for `cli`, all crates are published on `crates.io`.

## License

This project is licensed under the terms of the BSD 2-Clause license.
See `LICENSE.md` for more information.

[cli crate]: https://github.com/algbio/practical-omnitigs/tree/master/implementation/cli