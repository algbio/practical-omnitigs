use crate::omnitigs::{perform_linear_reduction, split_walks_at_node};
use crate::CliOptions;
use clap::Parser;
use compact_genome::implementation::DefaultSequenceStore;
use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use error_chain::bail;
use genome_graph::bigraph::interface::dynamic_bigraph::DynamicBigraph;
use genome_graph::bigraph::interface::static_bigraph::StaticEdgeCentricBigraph;
use genome_graph::bigraph::interface::BidirectedData;
use genome_graph::bigraph::traitgraph::interface::subgraph::SubgraphBase;
use genome_graph::types::{PetBCalm2EdgeGraph, PetWtdbg2DotGraph, PetWtdbg2Graph};
use log::info;
use omnitigs::omnitigs::remove_subwalks_and_reverse_complements_from_walks;
use omnitigs::omnitigs::Omnitigs;
use omnitigs::traitgraph::index::GraphIndex;
use omnitigs::traitgraph::interface::{GraphBase, ImmutableGraphContainer};
use omnitigs::traitgraph::walks::VecEdgeWalk;
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Write};
use traitsequence::interface::Sequence;

#[derive(Parser)]
pub struct ComputeMultiSafeWalksCommand {
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
        help = "The file the omnitigs are stored into in fasta format."
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
        long,
        help = "Connect all sources and sinks to a global node to make the graph strongly connected. Then compute omnitigs and split them if they contain the global node. If the graph does not become strongly connected, abort."
    )]
    pub linear_reduction: bool,

    #[clap(
        long,
        help = "When performing the linear reduction, report walks only on the SCC that contains the dummy node."
    )]
    pub linear_reduction_use_scc: bool,
}

fn print_multi_safe_walks_statistics<Graph: GraphBase>(
    maximal_multi_safe_walks: &Omnitigs<Graph>,
    latex_file: &mut Option<BufWriter<File>>,
) -> crate::Result<()> {
    info!("");
    info!(" === Multi-Safe Walk Statistics === ");
    info!("");

    if !maximal_multi_safe_walks.omnitigs_per_macrotig().is_empty() {
        let min_omnitigs_per_macrotig = maximal_multi_safe_walks
            .omnitigs_per_macrotig()
            .iter()
            .min()
            .unwrap();
        let max_omnitigs_per_macrotig = maximal_multi_safe_walks
            .omnitigs_per_macrotig()
            .iter()
            .max()
            .unwrap();
        let median_omnitigs_per_macrotig =
            statistical::median(maximal_multi_safe_walks.omnitigs_per_macrotig());
        let mean_omnitigs_per_macrotig = statistical::mean(
            &maximal_multi_safe_walks
                .omnitigs_per_macrotig()
                .iter()
                .map(|i| *i as f64)
                .collect::<Vec<_>>(),
        );

        info!(
            "Minimum non-trivial multi-safe walks per macrotig: {}",
            min_omnitigs_per_macrotig
        );
        info!(
            "Maximum non-trivial multi-safe walks per macrotig: {}",
            max_omnitigs_per_macrotig
        );
        info!(
            "Median non-trivial multi-safe walks per macrotig: {}",
            median_omnitigs_per_macrotig
        );
        info!(
            "Mean non-trivial multi-safe walks per macrotig: {:.1}",
            mean_omnitigs_per_macrotig
        );

        if let Some(latex_file) = latex_file.as_mut() {
            writeln!(
                latex_file,
                "min non-trivial multi-safe walks per macrotig & {} \\\\",
                min_omnitigs_per_macrotig
            )?;
            writeln!(
                latex_file,
                "max non-trivial multi-safe walks per macrotig & {} \\\\",
                max_omnitigs_per_macrotig
            )?;
            writeln!(
                latex_file,
                "median non-trivial multi-safe walks per macrotig & {} \\\\",
                median_omnitigs_per_macrotig
            )?;
            writeln!(
                latex_file,
                "mean non-trivial multi-safe walks per macrotig & {:.1} \\\\",
                mean_omnitigs_per_macrotig
            )?;
        }
    }

    let min_omnitig_len = maximal_multi_safe_walks
        .iter()
        .map(Sequence::len)
        .min()
        .unwrap();
    let max_omnitig_len = maximal_multi_safe_walks
        .iter()
        .map(Sequence::len)
        .max()
        .unwrap();
    let median_omnitigs_len = statistical::median(
        &maximal_multi_safe_walks
            .iter()
            .map(Sequence::len)
            .collect::<Vec<_>>(),
    );
    let mean_omnitig_len = statistical::mean(
        &maximal_multi_safe_walks
            .iter()
            .map(|o| o.len() as f64)
            .collect::<Vec<_>>(),
    );
    info!("Minimum edge length: {}", min_omnitig_len);
    info!("Maximum edge length: {}", max_omnitig_len);
    info!("Median edge length: {}", median_omnitigs_len);
    info!("Mean edge length: {:.1}", mean_omnitig_len);

    if let Some(latex_file) = latex_file.as_mut() {
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

fn compute_multi_safe_walks_from_graph<
    Graph: DynamicBigraph + Clone + Default + StaticEdgeCentricBigraph + SubgraphBase<RootGraph = Graph>,
>(
    genome_graph: &Graph,
    compute_multi_safe_walks_command: &ComputeMultiSafeWalksCommand,
    latex_file: &mut Option<BufWriter<File>>,
) -> crate::Result<Vec<VecEdgeWalk<Graph>>>
where
    Graph::NodeData: Default + Clone + Eq + Debug,
    Graph::EdgeData: Default + Clone + BidirectedData + Eq,
    <Graph as ImmutableGraphContainer>::NodeIndicesCopied: DoubleEndedIterator + ExactSizeIterator,
{
    debug_assert!(genome_graph.verify_edge_mirror_property());

    Ok(if compute_multi_safe_walks_command.linear_reduction {
        let (omnitig_graph, unreduce_map) = perform_linear_reduction(
            genome_graph,
            compute_multi_safe_walks_command.linear_reduction_use_scc,
        )?;

        info!("Computing maximal multi-safe walks");
        let mut omnitigs = Omnitigs::compute(&omnitig_graph);
        omnitigs.remove_reverse_complements(&omnitig_graph);
        print_multi_safe_walks_statistics(&omnitigs, latex_file)?;
        let mut omnitigs: Vec<_> = omnitigs.into_iter().map(Into::into).collect();
        split_walks_at_node(
            &mut omnitigs,
            omnitig_graph.node_indices().last().unwrap(),
            &omnitig_graph,
        )?;

        info!("Computing additional trivial multi-safe walks");
        let mut trivial_omnitigs = Omnitigs::compute_trivial_only_non_scc(genome_graph);
        trivial_omnitigs.remove_reverse_complements(genome_graph);

        info!("Merging tigs");
        for omnitig in &mut omnitigs {
            for edge in omnitig {
                *edge = unreduce_map[edge.as_usize()];
            }
        }
        omnitigs.extend(trivial_omnitigs.into_iter().map(Into::into));
        remove_subwalks_and_reverse_complements_from_walks(&mut omnitigs, genome_graph);
        omnitigs
    } else {
        info!("Computing maximal multi-safe walks");
        let mut omnitigs = Omnitigs::compute(genome_graph);
        omnitigs.remove_reverse_complements(genome_graph);
        print_multi_safe_walks_statistics(&omnitigs, latex_file)?;
        omnitigs.into_iter().map(Into::into).collect()
    })
}

pub(crate) fn compute_omnitigs(
    _options: &CliOptions,
    subcommand: &ComputeMultiSafeWalksCommand,
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
            info!(
                "Graph has {} nodes and {} edges",
                genome_graph.node_count(),
                genome_graph.edge_count()
            );

            let omnitigs =
                compute_multi_safe_walks_from_graph(&genome_graph, subcommand, &mut latex_file)?;

            info!(
                "Storing maximal multi-safe walks as fasta to '{}'",
                subcommand.output
            );
            genome_graph::io::fasta::write_walks_as_fasta_file(
                &genome_graph,
                &sequence_store,
                kmer_size,
                omnitigs.iter(),
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
            info!(
                "Graph has {} nodes and {} edges",
                genome_graph.node_count(),
                genome_graph.edge_count()
            );

            let omnitigs =
                compute_multi_safe_walks_from_graph(&genome_graph, subcommand, &mut latex_file)?;

            if subcommand.output_as_wtdbg2_node_ids {
                info!(
                    "Storing multi-safe walks as node ids to '{}'",
                    subcommand.output
                );
                genome_graph::io::wtdbg2::write_contigs_as_wtdbg2_node_ids_to_file(
                    &genome_graph,
                    omnitigs.iter(),
                    &subcommand.output,
                )?;
            } else {
                info!(
                    "Storing multi-safe walks as .ctg.lay to '{}'",
                    subcommand.output
                );
                genome_graph::io::wtdbg2::write_contigs_to_wtdbg2_to_file(
                    &genome_graph,
                    omnitigs.iter(),
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

            info!(
                "Graph has {} nodes and {} edges",
                genome_graph.node_count(),
                genome_graph.edge_count()
            );

            let omnitigs =
                compute_multi_safe_walks_from_graph(&genome_graph, subcommand, &mut latex_file)?;

            info!(
                "Storing multi-safe walks as node ids to '{}'",
                subcommand.output
            );
            genome_graph::io::wtdbg2::dot::write_dot_contigs_as_wtdbg2_node_ids_to_file(
                &genome_graph,
                omnitigs.iter(),
                &subcommand.output,
            )?;
        }

        unknown => bail!("Unknown file format: {}", unknown),
    }

    Ok(())
}
