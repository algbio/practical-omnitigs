#![recursion_limit = "1024"]
#[macro_use]
extern crate error_chain;

use clap::Clap;
use error_chain::{ChainedError, ExitCode};
use genome_graph::types::PetBCalm2Graph;

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

fn main() {
    use std::io::Write;

    ::std::process::exit(match run() {
        Ok(()) => ExitCode::code(()),
        Err(ref e) => {
            write!(
                &mut ::std::io::stderr(),
                "{}",
                ChainedError::display_chain(e)
            )
            .expect("Error writing to stderr");

            1
        }
    });
}

fn run() -> Result<()> {
    let options = &CliOptions::parse();

    println!("Hello");

    match &options.subcommand {
        Command::Verify => verify(options),
    }?;

    println!("Goodbye");
    Ok(())
}

fn verify(options: &CliOptions) -> Result<()> {
    println!(
        "Input file: {}\nOutput file: {}",
        options.input,
        options.output.as_ref().unwrap_or(&"None".to_owned())
    );

    let genome_graph: PetBCalm2Graph =
        genome_graph::io::bcalm2::load_bigraph_from_bcalm2(&options.input)?;
    println!("{:?}", genome_graph);
    Ok(())
}
