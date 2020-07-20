use clap::Clap;
use genome_graph::bigraph::*;

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
    let options = &CliOptions::parse();

    println!("Hello");

    match &options.subcommand {
        Command::Verify => verify(options),
    }

    println!("Goodbye");
}

fn verify(options: &CliOptions) {
    println!(
        "Input file: {}\nOutput file: {}",
        options.input,
        options.output.as_ref().unwrap_or(&"None".to_owned())
    );

    let genome_graph: NodeBigraphWrapper<genome_graph::BCalm2NodeData, (), usize, genome_graph::bigraph::petgraph::Graph<genome_graph::BCalm2NodeData, (), genome_graph::bigraph::petgraph::Directed, usize>> = genome_graph::load_bigraph_from_bcalm2(&options.input).unwrap();
    println!("{:?}", genome_graph);
}
