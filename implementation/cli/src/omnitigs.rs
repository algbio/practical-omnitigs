use crate::CliOptions;
use clap::Clap;
use genome_graph::types::PetBCalm2EdgeGraph;
use omnitigs::omnitigs::Omnitigs;
use std::io::Write;
use traitsequence::interface::Sequence;

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

    #[clap(
        short,
        long,
        about = "A file to output the properties and statistics computed by this command formatted as a LaTeX table"
    )]
    pub latex: Option<String>,
}

pub(crate) fn compute_omnitigs(
    options: &CliOptions,
    subcommand: &ComputeOmnitigsCommand,
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

    info!("Computing maximal omnitigs");
    let mut maximal_omnitigs = Omnitigs::compute(&genome_graph);
    info!("Removing reverse complements");
    maximal_omnitigs.remove_reverse_complements(&genome_graph);

    info!("");
    info!(" === Omnitig Statistics === ");
    info!("");

    let min_omnitigs_per_macrotig = maximal_omnitigs
        .omnitigs_per_macrotig()
        .iter()
        .min()
        .unwrap();
    let max_omnitigs_per_macrotig = maximal_omnitigs
        .omnitigs_per_macrotig()
        .iter()
        .max()
        .unwrap();
    let median_omnitigs_per_macrotig =
        statistical::median(maximal_omnitigs.omnitigs_per_macrotig());
    let mean_omnitigs_per_macrotig = statistical::mean(
        &maximal_omnitigs
            .omnitigs_per_macrotig()
            .iter()
            .map(|i| *i as f64)
            .collect::<Vec<_>>(),
    );
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

    info!(
        "Minimum non-trivial omnitigs per macrotig: {}",
        min_omnitigs_per_macrotig
    );
    info!(
        "Maximum non-trivial omnitigs per macrotig: {}",
        max_omnitigs_per_macrotig
    );
    info!(
        "Median non-trivial omnitigs per macrotig: {}",
        median_omnitigs_per_macrotig
    );
    info!(
        "Mean non-trivial omnitigs per macrotig: {}",
        mean_omnitigs_per_macrotig
    );
    info!("Minimum edge length: {}", min_omnitig_len);
    info!("Maximum edge length: {}", max_omnitig_len);
    info!("Median edge length: {}", median_omnitigs_len);
    info!("Mean edge length: {}", mean_omnitig_len);

    if let Some(latex_file) = &mut latex_file {
        writeln!(
            latex_file,
            "min non-trivial omnitigs per macrotig & {} \\\\",
            min_omnitigs_per_macrotig
        )?;
        writeln!(
            latex_file,
            "max non-trivial omnitigs per macrotig & {} \\\\",
            max_omnitigs_per_macrotig
        )?;
        writeln!(
            latex_file,
            "median non-trivial omnitigs per macrotig & {} \\\\",
            median_omnitigs_per_macrotig
        )?;
        writeln!(
            latex_file,
            "mean non-trivial omnitigs per macrotig & {} \\\\",
            mean_omnitigs_per_macrotig
        )?;
        writeln!(latex_file, "min edge length & {} \\\\", min_omnitig_len)?;
        writeln!(latex_file, "max edge length & {} \\\\", max_omnitig_len)?;
        writeln!(
            latex_file,
            "median edge length & {} \\\\",
            median_omnitigs_len
        )?;
        writeln!(latex_file, "mean edge length & {} \\\\", mean_omnitig_len)?;
    }

    info!("");
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
