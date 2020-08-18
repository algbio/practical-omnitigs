#![recursion_limit = "1024"]
#[macro_use]
extern crate error_chain;
#[macro_use]
extern crate log;

use clap::Clap;
use error_chain::{ChainedError, ExitCode};
use simplelog::{CombinedLogger, Config, LevelFilter, TermLogger, TerminalMode};

mod verify;

error_chain! {
    links {
        GenomeGraph(genome_graph::error::Error, genome_graph::error::ErrorKind);
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

    #[clap(
        short,
        long,
        about = "The output file, depending on the subcommand used"
    )]
    pub output: Option<String>,

    #[clap(subcommand)]
    pub subcommand: Command,
}

#[derive(Clap)]
enum Command {
    #[clap(
        about = "Prints statistics about the input graph, and saves it back to disc for verification purposes if --output is given"
    )]
    Verify,
    #[clap(
        about = "Same as verify, but loads the input graph node-centric instead of edge-centric"
    )]
    VerifyNodeCentric,
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
        Command::Verify => verify::verify_edge_centric(options),
        Command::VerifyNodeCentric => verify::verify_node_centric(options),
    }?;

    info!("Goodbye");
    Ok(())
}
