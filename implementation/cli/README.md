# Practical Omnitigs CLI

The command line interface of the omnitig related algorithms.

## Usage

In the folder that contains this README, run the following command.
```commandline
cargo run --release -- <...>
```

Replace `<...>` with the subcommands and arguments for the CLI.
For example:
```commandline
cargo run --release -- --input <your_bcalm_output_file.unitigs.fa> verify --kmer-size <your bcalm2 kmer size>
```

Use the following to get a help message with detailed usage information.
```commandline
cargo run --release -- --help
```