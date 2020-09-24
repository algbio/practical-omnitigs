#![recursion_limit = "1024"]
#[macro_use]
extern crate error_chain;
#[macro_use]
extern crate log;

use clap::Clap;
use error_chain::{ChainedError, ExitCode};
use simplelog::{CombinedLogger, Config, LevelFilter, TermLogger, TerminalMode};

mod circularise_records;
mod filter;
mod omnitigs;
mod trivial_omnitigs;
mod verify;
mod verify_genome;

error_chain! {
    foreign_links {
        Io(std::io::Error);
    }

    links {
        GenomeGraph(genome_graph::error::Error, genome_graph::error::ErrorKind);
        VerifyGenome(verify_genome::Error, verify_genome::ErrorKind);
    }

    errors {
        Parameter {
            description("a parameter was missing, superfluous or had an illegal value, see the log for more details")
            display("a parameter was missing, superfluous or had an illegal value, see the log for more details")
        }
    }
}

#[derive(Clap)]
#[clap(name = "Practical Omnitigs", version = env!("CARGO_PKG_VERSION"), author = "Sebastian Schmidt <sebastian.schmidt@helsinki.fi>")]
struct CliOptions {
    #[clap(
        short,
        long,
        about = "The input file, depending on the subcommand used"
    )]
    pub input: String,

    #[clap(subcommand)]
    pub subcommand: Command,
}

#[derive(Clap)]
enum Command {
    #[clap(
        about = "Prints statistics about the input graph, and saves it back to disc for verification purposes if --output is given."
    )]
    Verify(verify::VerifyEdgeCentricCommand),
    #[clap(
        about = "Same as verify, but loads the input graph node-centric instead of edge-centric."
    )]
    VerifyNodeCentric(verify::VerifyNodeCentricCommand),
    #[clap(about = "Verifies that no record of the genome contains illegal characters.")]
    VerifyGenome,
    #[clap(about = "Circularises each record in a genome by appending a prefix to itself.")]
    CirculariseGenome(circularise_records::CirculariseGenomeCommand),
    #[clap(about = "Filters records from a fasta file.")]
    Filter(filter::FilterCommand),
    #[clap(about = "Computes the maximal omnitigs of the input graph.")]
    ComputeOmnitigs(omnitigs::ComputeOmnitigsCommand),
    #[clap(about = "Computes the maximal trivial omnitigs of the input graph.")]
    ComputeTrivialOmnitigs(trivial_omnitigs::ComputeTrivialOmnitigsCommand),
}

// The main is unpacked from an error-chain macro.
// Using just the macro makes IntelliJ complain that there would be no main.
// The real main (programmed manually) is run, below this method.
fn main() {
    ::std::process::exit(match run() {
        Ok(()) => ExitCode::code(()),
        Err(ref e) => {
            error!("{}", ChainedError::display_chain(e));
            1
        }
    });
}

fn initialise_logging() {
    CombinedLogger::init(vec![TermLogger::new(
        LevelFilter::Trace,
        Config::default(),
        TerminalMode::Mixed,
    )])
    .unwrap();

    info!("Logging initialised successfully");
}

fn run() -> Result<()> {
    let options = &CliOptions::parse();
    initialise_logging();

    info!("Hello");

    match &options.subcommand {
        Command::Verify(subcommand) => verify::verify_edge_centric(options, subcommand),
        Command::VerifyNodeCentric(subcommand) => verify::verify_node_centric(options, subcommand),
        Command::VerifyGenome => verify_genome::verify_genome(options),
        Command::CirculariseGenome(subcommand) => {
            circularise_records::circularise_records(options, subcommand)
        }
        Command::Filter(subcommand) => filter::filter_records(options, subcommand),
        Command::ComputeOmnitigs(subcommand) => omnitigs::compute_omnitigs(options, subcommand),
        Command::ComputeTrivialOmnitigs(subcommand) => {
            trivial_omnitigs::compute_trivial_omnitigs(options, subcommand)
        }
    }?;

    info!("Goodbye");
    Ok(())
}
