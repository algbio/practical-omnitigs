use crate::CliOptions;
use clap::Clap;
use genome_graph::types::PetBCalm2EdgeGraph;
use omnitigs::omnitigs::Omnitigs;

#[derive(Clap)]
pub struct ComputeOmnitigsCommand {
    #[clap(
        short,
        long,
        about = "The kmer size selected when generating the input with bcalm2"
    )]
    pub kmer_size: usize,

    #[clap(
        short,
        long,
        about = "The file the omnitigs are stored into in fasta format."
    )]
    pub output: String,
}

pub(crate) fn compute_omnitigs(
    options: &CliOptions,
    subcommand: &ComputeOmnitigsCommand,
) -> crate::Result<()> {
    info!(
        "Reading bigraph from '{}' with kmer size {}",
        options.input, subcommand.kmer_size
    );
    let genome_graph: PetBCalm2EdgeGraph =
        genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file(
            &options.input,
            subcommand.kmer_size,
        )?;

    info!("Computing maximal omnitigs");
    let maximal_omnitigs = Omnitigs::compute(&genome_graph);
    info!(
        "Storing maximal omnitigs as fasta to '{}'",
        subcommand.output
    );
    genome_graph::io::fasta::write_walks_as_fasta_file(
        &genome_graph,
        subcommand.kmer_size,
        maximal_omnitigs.iter(),
        &subcommand.output,
    )?;

    Ok(())
}
