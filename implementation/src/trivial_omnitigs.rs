use crate::CliOptions;
use clap::Parser;
use compact_genome::implementation::{DefaultSequenceStore, DefaultSequenceStoreHandle};
use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use error_chain::{bail, ensure};
use genome_graph::io::gfa::PetGfaGraph;
use genome_graph::types::{PetBCalm2EdgeGraph, PetWtdbg2DotGraph, PetWtdbg2Graph};
use log::{info, warn};
use omnitigs::omnitigs::{NodeCentricOmnitigs, Omnitigs, TruncationMode};
use omnitigs::traitgraph::interface::GraphBase;
use omnitigs::traitgraph::interface::ImmutableGraphContainer;
use std::fs::File;
use std::io::{BufWriter, Write};
use traitgraph_algo::components::is_strongly_connected;
use traitsequence::interface::Sequence;

#[derive(Parser)]
pub struct ComputeTrivialOmnitigsCommand {
    #[clap(
        short,
        long,
        default_value = "bcalm2",
        help = "The format of the input and output files. If bcalm2, the input file is in bcalm2 format and the output file is in fasta format. If wtdbg2, the inputs are .1.nodes and the .1.reads file and the reads file from which these were generated, and the output is the .ctg.lay file. If dot, then the input is a .dot file and the output is a list of sequences of node ids. If hifiasm, the input is in gfa format, and the output in fasta format."
    )]
    pub file_format: String,

    #[clap(short, long, help = "The input files in the specified format")]
    pub input: Vec<String>,

    #[clap(
        short,
        long,
        help = "The kmer size selected when generating the input with bcalm2"
    )]
    pub kmer_size: Option<usize>,

    #[clap(
        short,
        long,
        help = "The file the trivial omnitigs are stored into in fasta format"
    )]
    pub output: String,

    #[clap(
        short,
        long,
        help = "A file to output the properties and statistics computed by this command formatted as a LaTeX table"
    )]
    pub latex: Option<String>,

    #[clap(
        long,
        help = "Instead of outputting unitigs as .ctg.lay file, output them as sequences of node ids"
    )]
    pub output_as_wtdbg2_node_ids: bool,

    #[clap(
        short,
        long,
        help = "Use algorithms that handle not strongly connected graphs, but are slower"
    )]
    pub non_scc: bool,

    #[clap(
        short,
        long,
        help = "Truncate univocal extensions that are not multi-safe because they repeat edges"
    )]
    pub multi_safe: bool,

    #[clap(
        long,
        help = "Do not remove reverse complements from the resulting trivial omnitigs"
    )]
    pub keep_reverse_complements: bool,
}

fn print_trivial_omnitigs_statistics<Graph: GraphBase>(
    maximal_omnitigs: &Omnitigs<Graph>,
    latex_file: &mut Option<BufWriter<File>>,
) -> crate::Result<()> {
    info!("");
    info!(" === Trivial Omnitig Statistics === ");
    info!("");

    if maximal_omnitigs.is_empty() {
        warn!("Found no trivial omnitigs");
    }

    let omnitig_lengths_usize = maximal_omnitigs
        .iter()
        .map(Sequence::len)
        .collect::<Vec<_>>();
    let omnitig_lengths_f64 = maximal_omnitigs
        .iter()
        .map(|o| o.len() as f64)
        .collect::<Vec<_>>();

    let min_omnitig_len = maximal_omnitigs
        .iter()
        .map(Sequence::len)
        .min()
        .unwrap_or(0);
    let max_omnitig_len = maximal_omnitigs
        .iter()
        .map(Sequence::len)
        .max()
        .unwrap_or(0);
    let median_omnitigs_len = if omnitig_lengths_usize.is_empty() {
        0
    } else {
        statistical::median(&omnitig_lengths_usize)
    };
    let mean_omnitig_len = if omnitig_lengths_f64.is_empty() {
        0.0
    } else {
        statistical::mean(&omnitig_lengths_f64)
    };

    info!("Minimum edge length: {}", min_omnitig_len);
    info!("Maximum edge length: {}", max_omnitig_len);
    info!("Median edge length: {}", median_omnitigs_len);
    info!("Mean edge length: {:.1}", mean_omnitig_len);

    if let Some(latex_file) = latex_file.as_mut() {
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
    Ok(())
}

pub(crate) fn compute_trivial_omnitigs(
    _options: &CliOptions,
    subcommand: &ComputeTrivialOmnitigsCommand,
) -> crate::Result<()> {
    let mut latex_file = if let Some(latex_file_name) = &subcommand.latex {
        info!("Creating/truncating LaTeX file '{}'", latex_file_name);
        Some(std::io::BufWriter::new(std::fs::File::create(
            latex_file_name,
        )?))
    } else {
        None
    };

    let mut sequence_store = DefaultSequenceStore::default();

    match subcommand.file_format.as_str() {
        "bcalm2" => {
            if subcommand.output_as_wtdbg2_node_ids {
                bail!("Output as wtdbg2 node ids not supported for bcalm2 format");
            }

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
            let genome_graph: PetBCalm2EdgeGraph<_> =
                genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file(
                    input,
                    &mut sequence_store,
                    kmer_size,
                )?;

            info!("Computing maximal trivial omnitigs");
            let mut trivial_omnitigs = if subcommand.non_scc {
                Omnitigs::compute_trivial_only_non_scc(&genome_graph)
            } else {
                ensure!(is_strongly_connected(&genome_graph), "The graph is not strongly connected, but algorithms for not strongly connected graphs were not selected. Use --non-scc.");
                Omnitigs::compute_trivial_only(&genome_graph)
            };
            if !subcommand.keep_reverse_complements {
                info!("Removing reverse complements");
                trivial_omnitigs.remove_reverse_complements(&genome_graph);
            }
            assert!(!trivial_omnitigs.is_empty(), "Found no trivial omnitigs");

            if subcommand.multi_safe {
                trivial_omnitigs.transform_to_multi_safe_strict_model(TruncationMode::All);
            }

            print_trivial_omnitigs_statistics(&trivial_omnitigs, &mut latex_file)?;

            info!(
                "Storing maximal trivial omnitigs as fasta to '{}'",
                subcommand.output
            );
            genome_graph::io::fasta::write_walks_as_fasta_file(
                &genome_graph,
                &sequence_store,
                kmer_size,
                trivial_omnitigs.iter(),
                &subcommand.output,
            )?;
        }

        "hifiasm" => {
            if subcommand.output_as_wtdbg2_node_ids {
                bail!("Output as wtdbg2 node ids not supported for hifiasm format");
            }

            if subcommand.multi_safe {
                bail!("Input format hifiasm does not support multi-safe walks");
            }

            let input = if let Some(input) = subcommand.input.first() {
                input
            } else {
                bail!("No input file given")
            };
            info!("Reading bigraph from '{}'", input);
            let (genome_graph, _): (
                PetGfaGraph<(), (), DefaultSequenceStoreHandle<DnaAlphabet>>,
                _,
            ) = genome_graph::io::gfa::read_gfa_as_bigraph_from_file(
                input,
                &mut sequence_store,
                true,
                true,
            )?;

            info!(
                "Graph has {} nodes and {} edges",
                genome_graph.node_count(),
                genome_graph.edge_count()
            );

            info!("Computing maximal trivial omnitigs");
            let mut trivial_omnitigs = if subcommand.non_scc {
                Vec::compute_trivial_node_centric_omnitigs_non_scc(&genome_graph)
            } else {
                ensure!(is_strongly_connected(&genome_graph), "The graph is not strongly connected, but algorithms for not strongly connected graphs were not selected. Use --non-scc.");
                Vec::compute_trivial_node_centric_omnitigs(&genome_graph)
            };
            if !subcommand.keep_reverse_complements {
                info!("Removing reverse complements");
                trivial_omnitigs.remove_reverse_complements(&genome_graph);
            }

            //print_trivial_node_centric_omnitigs_statistics(&maximal_omnitigs, &mut latex_file)?;

            info!(
                "Storing maximal trivial omnitigs as fasta to '{}'",
                subcommand.output
            );
            genome_graph::io::fasta::write_node_centric_walks_with_variable_overlaps_as_fasta_file(
                &genome_graph,
                &sequence_store,
                trivial_omnitigs.iter(),
                &subcommand.output,
            )?;
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

            info!("Computing maximal trivial omnitigs");
            let mut trivial_omnitigs = if subcommand.non_scc {
                Omnitigs::compute_trivial_only_non_scc(&genome_graph)
            } else {
                ensure!(is_strongly_connected(&genome_graph), "The graph is not strongly connected, but algorithms for not strongly connected graphs were not selected. Use --non-scc.");
                Omnitigs::compute_trivial_only(&genome_graph)
            };
            if !subcommand.keep_reverse_complements {
                info!("Removing reverse complements");
                trivial_omnitigs.remove_reverse_complements(&genome_graph);
            }

            if subcommand.multi_safe {
                trivial_omnitigs.transform_to_multi_safe_strict_model(TruncationMode::All);
            }

            print_trivial_omnitigs_statistics(&trivial_omnitigs, &mut latex_file)?;

            if subcommand.output_as_wtdbg2_node_ids {
                info!(
                    "Storing trivial omnitigs as node ids to '{}'",
                    subcommand.output
                );
                genome_graph::io::wtdbg2::write_contigs_as_wtdbg2_node_ids_to_file(
                    &genome_graph,
                    trivial_omnitigs.iter(),
                    &subcommand.output,
                )?;
            } else {
                info!(
                    "Storing trivial omnitigs as .ctg.lay to '{}'",
                    subcommand.output
                );
                genome_graph::io::wtdbg2::write_contigs_to_wtdbg2_to_file(
                    &genome_graph,
                    trivial_omnitigs.iter(),
                    raw_reads_file,
                    &subcommand.output,
                )?;
            }
        }

        "dot" => {
            let dot_file = if let Some(file) = subcommand.input.iter().find(|f| f.ends_with(".dot"))
            {
                file
            } else {
                bail!("Missing .dot file")
            };
            info!("Reading bigraph from '{}'", dot_file);

            let genome_graph: PetWtdbg2DotGraph =
                genome_graph::io::wtdbg2::dot::read_graph_from_wtdbg2_dot_from_file(dot_file)?;

            info!("Computing maximal trivial omnitigs");
            let mut trivial_omnitigs = if subcommand.non_scc {
                Omnitigs::compute_trivial_only_non_scc(&genome_graph)
            } else {
                ensure!(is_strongly_connected(&genome_graph), "The graph is not strongly connected, but algorithms for not strongly connected graphs were not selected. Use --non-scc.");
                Omnitigs::compute_trivial_only(&genome_graph)
            };
            if !subcommand.keep_reverse_complements {
                info!("Removing reverse complements");
                trivial_omnitigs.remove_reverse_complements(&genome_graph);
            }

            if subcommand.multi_safe {
                trivial_omnitigs.transform_to_multi_safe_strict_model(TruncationMode::All);
            }

            print_trivial_omnitigs_statistics(&trivial_omnitigs, &mut latex_file)?;

            info!(
                "Storing trivial omnitigs as node ids to '{}'",
                subcommand.output
            );
            genome_graph::io::wtdbg2::dot::write_dot_contigs_as_wtdbg2_node_ids_to_file(
                &genome_graph,
                trivial_omnitigs.iter(),
                &subcommand.output,
            )?;
        }

        unknown => bail!("Unknown file format: {}", unknown),
    }

    Ok(())
}
