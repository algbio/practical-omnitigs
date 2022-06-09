use crate::CliOptions;
use clap::Parser;
use compact_genome::implementation::DefaultSequenceStore;
use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use genome_graph::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
use genome_graph::bigraph::traitgraph::interface::GraphBase;
use genome_graph::io::wtdbg2::build_wtdbg2_unitigs_graph;
use genome_graph::types::{PetBCalm2EdgeGraph, PetWtdbg2DotGraph, PetWtdbg2Graph};
use omnitigs::traitgraph::implementation::petgraph_impl::PetGraph;
use omnitigs::unitigs::EdgeUnitigs;
use std::fs::File;
use std::io::{BufWriter, Write};
use traitgraph_algo::components::{
    decompose_strongly_connected_components, decompose_weakly_connected_components,
};
use traitsequence::interface::Sequence;

#[derive(Parser)]
pub struct ComputeUnitigsCommand {
    #[clap(
        short,
        long,
        default_value = "bcalm2",
        help = "The format of the input and output files. If bcalm2, the input file is in bcalm2 format and the output file is in fasta format. If wtdbg2, the inputs are .1.nodes and the .1.reads file and the reads file from which these were generated, and the output is the .ctg.lay file. If dot, then the input is a .dot file and the output is a list of sequences of node ids."
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
        help = "The file the unitigs are stored into in the specified format"
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
        help = "Compare the unitigs produced by our algorithm to wtdbg2's contigs"
    )]
    pub compare_with_wtdbg2_contigs: bool,
}

fn print_unitig_statistics<Graph: GraphBase>(
    unitigs: &EdgeUnitigs<Graph>,
    latex_file: &mut Option<BufWriter<File>>,
) -> crate::Result<()> {
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
    info!("Mean edge length: {:.1}", mean_unitig_len);

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
        writeln!(latex_file, "min edge length & {} \\\\", min_unitig_len)?;
        writeln!(latex_file, "max edge length & {} \\\\", max_unitig_len)?;
        writeln!(
            latex_file,
            "median edge length & {} \\\\",
            median_unitig_len
        )?;
        writeln!(latex_file, "mean edge length & {:.1} \\\\", mean_unitig_len)?;
    }

    info!("");
    Ok(())
}

pub(crate) fn compute_unitigs(
    _options: &CliOptions,
    subcommand: &ComputeUnitigsCommand,
) -> crate::Result<()> {
    let mut latex_file = if let Some(latex_file_name) = &subcommand.latex {
        info!("Creating/truncating LaTeX file '{}'", latex_file_name);
        Some(std::io::BufWriter::new(std::fs::File::create(
            latex_file_name,
        )?))
    } else {
        None
    };

    let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::default();

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

            info!("Computing maximal unitigs");
            let mut unitigs = EdgeUnitigs::compute(&genome_graph);
            info!("Removing reverse complements");
            unitigs.remove_reverse_complements(&genome_graph);

            print_unitig_statistics(&unitigs, &mut latex_file)?;

            info!("Storing unitigs as fasta to '{}'", subcommand.output);
            genome_graph::io::fasta::write_walks_as_fasta_file(
                &genome_graph,
                &sequence_store,
                kmer_size,
                unitigs.iter(),
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

            info!("Computing maximal unitigs");
            let mut unitigs = EdgeUnitigs::compute(&genome_graph);
            info!("Removing reverse complements");
            unitigs.remove_reverse_complements(&genome_graph);

            print_unitig_statistics(&unitigs, &mut latex_file)?;

            if subcommand.compare_with_wtdbg2_contigs {
                info!("Investigating differences between wtdbg2 and our unitigs");
                let wtdbg2_unitigs_file =
                    dot_file[..dot_file.len() - 6].to_owned() + ".wtdbg2.ctg.lay";
                info!("Loading wtdbg2 unitigs from '{}'", wtdbg2_unitigs_file);
                let mut wtdbg2_unitigs =
                    genome_graph::io::wtdbg2::read_wtdbg2_contigs_from_file(&wtdbg2_unitigs_file)?;
                info!("Converting our unitigs to .ctg.lay format");
                let mut our_unitigs =
                    genome_graph::io::wtdbg2::convert_walks_to_wtdbg2_contigs_with_file(
                        &genome_graph,
                        unitigs.iter(),
                        &raw_reads_file,
                    )?;
                info!("Sorting unitigs");
                wtdbg2_unitigs.sort_contigs_topologically();
                our_unitigs.sort_contigs_topologically();
                wtdbg2_unitigs.update_indices();
                our_unitigs.update_indices();
                info!(" =========================");
                info!(" === Comparing unitigs ===");
                info!(" =========================");
                wtdbg2_unitigs.compare_contigs(&our_unitigs);
                drop(our_unitigs);

                info!(" ==============================");
                info!(" === Analysing unitig graph ===");
                info!(" ==============================");
                let unitig_graph: NodeBigraphWrapper<PetGraph<_, _>> =
                    build_wtdbg2_unitigs_graph(&wtdbg2_unitigs);
                drop(wtdbg2_unitigs);

                let wccs = decompose_weakly_connected_components(&unitig_graph);
                info!("Unitig graph has {} wccs", wccs.len());
                let sccs = decompose_strongly_connected_components(&unitig_graph);
                info!("Unitig graph has {} sccs", sccs.len());
            }

            if subcommand.output_as_wtdbg2_node_ids {
                info!("Storing unitigs as node ids to '{}'", subcommand.output);
                genome_graph::io::wtdbg2::write_contigs_as_wtdbg2_node_ids_to_file(
                    &genome_graph,
                    unitigs.iter(),
                    &subcommand.output,
                )?;
            } else {
                info!("Storing unitigs as .ctg.lay to '{}'", subcommand.output);
                genome_graph::io::wtdbg2::write_contigs_to_wtdbg2_to_file(
                    &genome_graph,
                    unitigs.iter(),
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

            info!("Computing maximal unitigs");
            let mut unitigs = EdgeUnitigs::compute(&genome_graph);
            info!("Removing reverse complements");
            unitigs.remove_reverse_complements(&genome_graph);

            print_unitig_statistics(&unitigs, &mut latex_file)?;

            info!("Storing unitigs as node ids to '{}'", subcommand.output);
            genome_graph::io::wtdbg2::dot::write_dot_contigs_as_wtdbg2_node_ids_to_file(
                &genome_graph,
                unitigs.iter(),
                &subcommand.output,
            )?;
        }

        unknown => bail!("Unknown file format: {}", unknown),
    }

    Ok(())
}
