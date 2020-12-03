use crate::CliOptions;
use clap::Clap;
use genome_graph::types::{PetBCalm2EdgeGraph, PetWtdbg2Graph};
use omnitigs::omnitigs::Omnitigs;
use omnitigs::traitgraph::interface::GraphBase;
use std::fs::File;
use std::io::{BufWriter, Write};
use traitsequence::interface::Sequence;

#[derive(Clap)]
pub struct ComputeOmnitigsCommand {
    #[clap(
        short,
        long,
        default_value = "bcalm2",
        about = "The format of the input and output files. If bcalm2, the input file is in bcalm2 format and the output file is in fasta format. If wtdbg2, the inputs are .1.nodes and the .1.reads file and the reads file from which these were generated, and the output is the .ctg.lay file."
    )]
    pub file_format: String,

    #[clap(short, long, about = "The input files in the specified format")]
    pub input: Vec<String>,

    #[clap(
        short,
        long,
        about = "The kmer size selected when generating the input with bcalm2"
    )]
    pub kmer_size: Option<usize>,

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

fn print_omnitigs_statistics<Graph: GraphBase>(
    maximal_omnitigs: &Omnitigs<Graph>,
    latex_file: &mut Option<BufWriter<File>>,
) -> crate::Result<()> {
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
        "Mean non-trivial omnitigs per macrotig: {:.1}",
        mean_omnitigs_per_macrotig
    );
    info!("Minimum edge length: {}", min_omnitig_len);
    info!("Maximum edge length: {}", max_omnitig_len);
    info!("Median edge length: {}", median_omnitigs_len);
    info!("Mean edge length: {:.1}", mean_omnitig_len);

    if let Some(latex_file) = latex_file.as_mut() {
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
            "mean non-trivial omnitigs per macrotig & {:.1} \\\\",
            mean_omnitigs_per_macrotig
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
    Ok(())
}

pub(crate) fn compute_omnitigs(
    _options: &CliOptions,
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

    match subcommand.file_format.as_str() {
        "bcalm2" => {
            let input = if let Some(input) = subcommand.input.first() {
                input
            } else {
                bail!("No input file given")
            };
            let kmer_size = if let Some(kmer_size) = subcommand.kmer_size {
                kmer_size
            } else {
                bail!("No kmer size given")
            };
            info!(
                "Reading bigraph from '{}' with kmer size {}",
                input, kmer_size
            );
            let genome_graph: PetBCalm2EdgeGraph =
                genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file(
                    &input, kmer_size,
                )?;

            info!("Computing maximal omnitigs");
            let mut maximal_omnitigs = Omnitigs::compute(&genome_graph);
            info!("Removing reverse complements");
            maximal_omnitigs.remove_reverse_complements(&genome_graph);

            print_omnitigs_statistics(&maximal_omnitigs, &mut latex_file)?;

            info!(
                "Storing maximal omnitigs as fasta to '{}'",
                subcommand.output
            );
            genome_graph::io::fasta::write_walks_as_fasta_file(
                &genome_graph,
                kmer_size,
                maximal_omnitigs.iter(),
                &subcommand.output,
            )?;

            Ok(())
        }

        "wtdbg2" => {
            let nodes_file =
                if let Some(file) = subcommand.input.iter().find(|f| f.ends_with(".3.nodes")) {
                    file
                } else {
                    bail!("Missing .3.nodes file")
                };
            let reads_file =
                if let Some(file) = subcommand.input.iter().find(|f| f.ends_with(".3.reads")) {
                    file
                } else {
                    bail!("Missing .3.reads file")
                };
            let dot_file = if let Some(file) = subcommand.input.iter().find(|f| f.ends_with(".dot"))
            {
                file
            } else {
                bail!("Missing .dot file")
            };
            let raw_reads_file =
                if let Some(file) = subcommand.input.iter().find(|f| f.ends_with(".fa")) {
                    file
                } else {
                    bail!("Missing raw reads file ending on .fa")
                };
            info!(
                "Reading bigraph from '{}', '{}', '{}' and '{}'",
                nodes_file, reads_file, dot_file, raw_reads_file
            );

            let genome_graph: PetWtdbg2Graph =
                genome_graph::io::wtdbg2::read_graph_from_wtdbg2_from_files(
                    nodes_file, reads_file, dot_file,
                )?;

            info!("Computing maximal unitigs");
            let mut omnitigs = Omnitigs::compute(&genome_graph);
            info!("Removing reverse complements");
            omnitigs.remove_reverse_complements(&genome_graph);

            print_omnitigs_statistics(&omnitigs, &mut latex_file)?;

            info!("Storing unitigs as .ctg.lay to '{}'", subcommand.output);
            genome_graph::io::wtdbg2::write_contigs_to_wtdbg2_to_file(
                &genome_graph,
                omnitigs.iter(),
                raw_reads_file,
                &subcommand.output,
            )?;

            Ok(())
        }

        unknown => bail!("Unknown file format: {}", unknown),
    }
}
