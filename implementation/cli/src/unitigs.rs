use crate::CliOptions;
use clap::Clap;
use genome_graph::types::PetBCalm2EdgeGraph;
use omnitigs::unitigs::EdgeUnitigs;
use std::io::Write;
use traitsequence::interface::Sequence;

#[derive(Clap)]
pub struct ComputeUnitigsCommand {
    #[clap(
        short,
        long,
        about = "The kmer size selected when generating the input with bcalm2"
    )]
    pub kmer_size: usize,

    #[clap(
        short,
        long,
        about = "The file the unitigs are stored into in fasta format."
    )]
    pub output: String,

    #[clap(
        short,
        long,
        about = "A file to output the properties and statistics computed by this command formatted as a LaTeX table"
    )]
    pub latex: Option<String>,
}

pub(crate) fn compute_unitigs(
    options: &CliOptions,
    subcommand: &ComputeUnitigsCommand,
) -> crate::Result<()> {
    let mut latex_file = if let Some(latex_file_name) = &subcommand.latex {
        info!("Creating/truncating LaTeX file");
        Some(std::io::BufWriter::new(std::fs::File::create(
            latex_file_name,
        )?))
    } else {
        None
    };

    info!(
        "Reading bigraph from '{}' with kmer size {}",
        options.input, subcommand.kmer_size
    );
    let genome_graph: PetBCalm2EdgeGraph =
        genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file(
            &options.input,
            subcommand.kmer_size,
        )?;

    info!("Computing maximal unitigs");
    let mut unitigs = EdgeUnitigs::compute(&genome_graph);
    info!("Removing reverse complements");
    unitigs.remove_reverse_complements(&genome_graph);

    info!("");
    info!(" === Unitig Statistics === ");
    info!("");

    let min_unitig_len = unitigs.iter().map(Sequence::len).min().unwrap();
    let max_unitig_len = unitigs.iter().map(Sequence::len).max().unwrap();
    let median_unitig_len =
        statistical::median(&unitigs.iter().map(Sequence::len).collect::<Vec<_>>());
    let mean_unitig_len =
        statistical::mean(&unitigs.iter().map(|o| o.len() as f64).collect::<Vec<_>>());

    info!("Minimum edge length: {}", min_unitig_len);
    info!("Maximum edge length: {}", max_unitig_len);
    info!("Median edge length: {}", median_unitig_len);
    info!("Mean edge length: {}", mean_unitig_len);

    if let Some(latex_file) = &mut latex_file {
        writeln!(
            latex_file,
            "min non-trivial omnitigs per macrotig & N/A \\\\"
        )?;
        writeln!(
            latex_file,
            "max non-trivial omnitigs per macrotig & N/A \\\\"
        )?;
        writeln!(
            latex_file,
            "median non-trivial omnitigs per macrotig & N/A \\\\"
        )?;
        writeln!(
            latex_file,
            "mean non-trivial omnitigs per macrotig & N/A \\\\"
        )?;
        writeln!(latex_file, "min edge length & {} \\\\", min_unitig_len)?;
        writeln!(latex_file, "max edge length & {} \\\\", max_unitig_len)?;
        writeln!(
            latex_file,
            "median edge length & {} \\\\",
            median_unitig_len
        )?;
        writeln!(latex_file, "mean edge length & {} \\\\", mean_unitig_len)?;
    }

    info!("");
    info!("Storing unitigs as fasta to '{}'", subcommand.output);
    genome_graph::io::fasta::write_walks_as_fasta_file(
        &genome_graph,
        subcommand.kmer_size,
        unitigs.iter(),
        &subcommand.output,
    )?;

    Ok(())
}
