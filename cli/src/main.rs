use clap::Clap;

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
}
