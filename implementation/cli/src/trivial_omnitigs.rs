use crate::CliOptions;
use clap::Clap;
use genome_graph::types::PetBCalm2EdgeGraph;
use omnitigs::omnitigs::Omnitigs;
use std::io::Write;
use traitsequence::interface::Sequence;

#[derive(Clap)]
pub struct ComputeTrivialOmnitigsCommand {
    #[clap(short, long, about = "The input file in bcalm2 format")]
    pub input: String,

    #[clap(
        short,
        long,
        about = "The kmer size selected when generating the input with bcalm2"
    )]
    pub kmer_size: usize,

    #[clap(
        short,
        long,
        about = "The file the trivial omnitigs are stored into in fasta format."
    )]
    pub output: String,

    #[clap(
        short,
        long,
        about = "A file to output the properties and statistics computed by this command formatted as a LaTeX table"
    )]
    pub latex: Option<String>,
}

pub(crate) fn compute_trivial_omnitigs(
    _options: &CliOptions,
    subcommand: &ComputeTrivialOmnitigsCommand,
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
        subcommand.input, subcommand.kmer_size
    );
    let genome_graph: PetBCalm2EdgeGraph =
        genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file(
            &subcommand.input,
            subcommand.kmer_size,
        )?;

    info!("Computing maximal trivial omnitigs");
    let mut maximal_omnitigs = Omnitigs::compute_trivial_only(&genome_graph);
    info!("Removing reverse complements");
    maximal_omnitigs.remove_reverse_complements(&genome_graph);

    info!("");
    info!(" === Trivial Omnitig Statistics === ");
    info!("");

    let min_omnitig_len = maximal_omnitigs.iter().map(Sequence::len).min().unwrap();
    let max_omnitig_len = maximal_omnitigs.iter().map(Sequence::len).max().unwrap();
    let median_omnitigs_len = statistical::median(
        &maximal_omnitigs
            .iter()
            .map(Sequence::len)
            .collect::<Vec<_>>(),
    );
    let mean_omnitig_len = statistical::mean(
        &maximal_omnitigs
            .iter()
            .map(|o| o.len() as f64)
            .collect::<Vec<_>>(),
    );

    info!("Minimum edge length: {}", min_omnitig_len);
    info!("Maximum edge length: {}", max_omnitig_len);
    info!("Median edge length: {}", median_omnitigs_len);
    info!("Mean edge length: {:.1}", mean_omnitig_len);

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
        writeln!(latex_file, "min edge length & {} \\\\", min_omnitig_len)?;
        writeln!(latex_file, "max edge length & {} \\\\", max_omnitig_len)?;
        writeln!(
            latex_file,
            "median edge length & {} \\\\",
            median_omnitigs_len
        )?;
        writeln!(
            latex_file,
            "mean edge length & {:.1} \\\\",
            mean_omnitig_len
        )?;
    }

    info!("");
    info!(
        "Storing maximal trivial omnitigs as fasta to '{}'",
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
