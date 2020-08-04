#![recursion_limit = "1024"]
#[macro_use]
extern crate error_chain;
#[macro_use]
extern crate log;

use clap::Clap;
use colored::*;
use error_chain::{ChainedError, ExitCode};
use genome_graph::bigraph::ImmutableGraphContainer;
use genome_graph::types::PetBCalm2Graph;
use simplelog::{CombinedLogger, Config, LevelFilter, TermLogger, TerminalMode};

error_chain! {
    links {
        GenomeGraph(genome_graph::Error, genome_graph::ErrorKind);
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
    Verify,
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
        Command::Verify => verify(options),
    }?;

    info!("Goodbye");
    Ok(())
}

fn verify(options: &CliOptions) -> Result<()> {
    info!("Input file: {}", options.input);
    info!(
        "Output file: {}",
        options.output.as_ref().unwrap_or(&"None".to_owned())
    );

    let genome_graph: PetBCalm2Graph =
        genome_graph::io::bcalm2::read_bigraph_from_bcalm2_file(&options.input)?;

    info!("");
    info!("========================");
    info!("=== Graph Statistics ===");
    info!("========================");
    info!("");
    info!(
        "{} binodes, {} biedges",
        genome_graph.node_count() / 2,
        genome_graph.edge_count() / 2
    );
    if genome_graph.node_count() % 2 != 0 || genome_graph.edge_count() % 2 != 0 {
        error!("Internal error: uneven amount of nodes or edges in a bidirected graph");
    }

    let uncompacted_unitig_amount =
        genome_graph::bigraph::unitigs::count_uncompacted_unitigs(&genome_graph);
    if uncompacted_unitig_amount % 2 != 0 {
        error!("Internal error: uneven amount of uncompacted unitigs in a bidirected graph");
    }

    let log_string = format!(
        "{} uncompacted bidirected unitigs",
        uncompacted_unitig_amount / 2
    );
    if uncompacted_unitig_amount == 0 {
        info!("{}", log_string);
    } else {
        warn!("{}", log_string.yellow());
    }

    info!("");

    if let Some(output) = &options.output {
        genome_graph::io::bcalm2::write_bigraph_to_bcalm2_file(&genome_graph, output)?;
    }
    Ok(())
}
