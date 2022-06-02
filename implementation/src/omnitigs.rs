use crate::{CliOptions, ErrorKind};
use clap::Parser;
use compact_genome::implementation::DefaultSequenceStore;
use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use genome_graph::bigraph::interface::dynamic_bigraph::DynamicBigraph;
use genome_graph::types::{PetBCalm2EdgeGraph, PetWtdbg2DotGraph, PetWtdbg2Graph};
use omnitigs::omnitigs::remove_subwalks_and_reverse_complements_from_walks;
use omnitigs::omnitigs::Omnitigs;
use omnitigs::traitgraph::interface::{GraphBase, ImmutableGraphContainer, StaticGraph};
use omnitigs::traitgraph::walks::VecEdgeWalk;
use std::fs::File;
use std::io::{BufWriter, Write};
use traitgraph_algo::components::{decompose_strongly_connected_components, is_strongly_connected};
use traitsequence::interface::Sequence;

#[derive(Parser)]
pub struct ComputeOmnitigsCommand {
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
}

fn print_omnitigs_statistics<Graph: GraphBase>(
    maximal_omnitigs: &Omnitigs<Graph>,
    latex_file: &mut Option<BufWriter<File>>,
) -> crate::Result<()> {
    info!("");
    info!(" === Omnitig Statistics === ");
    info!("");

    if !maximal_omnitigs.omnitigs_per_macrotig().is_empty() {
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
        }
    }

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

fn perform_linear_reduction<Graph: DynamicBigraph>(graph: &mut Graph) -> crate::Result<()>
where
    Graph::NodeData: Default,
    Graph::EdgeData: Default,
{
    info!("Performing linear reduction");
    let global_node = graph.add_node(Default::default());
    graph.set_mirror_nodes(global_node, global_node);
    let mut inserted_edge_count = 0usize;

    for node in graph.node_indices().rev().skip(1).rev() {
        if graph.out_degree(node) == 0 {
            graph.add_edge(node, global_node, Default::default());
            inserted_edge_count += 1;
        }
        if graph.in_degree(node) == 0 {
            graph.add_edge(global_node, node, Default::default());
            inserted_edge_count += 1;
        }
    }

    if inserted_edge_count == 0 {
        graph.remove_node(global_node);
    }

    info!("Linear reduction inserted {} edges", inserted_edge_count);

    if !is_strongly_connected(graph) {
        error!("Graph is not strongly connected after linear reduction");
        return Err(ErrorKind::LinearReductionNotStronglyConnected.into());
    }

    Ok(())
}

fn split_walks_at_node<Graph: StaticGraph>(
    walks: &mut Vec<VecEdgeWalk<Graph>>,
    split_node: Graph::NodeIndex,
    graph: &Graph,
) -> crate::Result<()> {
    info!("Splitting walks to revert linear reduction");
    // also check the newly added walks, because they might contain another instance of node
    let mut index = 0;
    while index < walks.len() {
        let walk = &mut walks[index];
        let mut removed_walk = false;
        while graph.edge_endpoints(*walk.first().unwrap()).from_node == split_node && !removed_walk
        {
            walk.remove(0);
            if walk.is_empty() {
                walks.remove(index);
                removed_walk = true;
                break;
            }
        }
        let walk = &mut walks[index];
        while graph.edge_endpoints(*walk.last().unwrap()).to_node == split_node && !removed_walk {
            walk.pop();
            if walk.is_empty() {
                walks.remove(index);
                removed_walk = true;
                break;
            }
        }
        if removed_walk {
            continue;
        }
        let walk = &mut walks[index];

        // here, only inner nodes can be the split node
        if let Some((walk_index, _)) = walk
            .iter()
            .enumerate()
            .skip(1)
            .find(|(_, &edge_index)| graph.edge_endpoints(edge_index).from_node == split_node)
        {
            debug!("Splitting walk {} at {}", index, walk_index);
            // arcs incident to split nodes are part of the linear reduction as well
            let right_half = walk.split_off(walk_index + 1);
            walk.pop();
            walk.pop();
            if !right_half.is_empty() {
                walks.push(right_half);
            }
        }

        index += 1;
    }

    Ok(())
}

pub(crate) fn compute_omnitigs(
    _options: &CliOptions,
    subcommand: &ComputeOmnitigsCommand,
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
            let mut genome_graph: PetBCalm2EdgeGraph<_> =
                genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file(
                    &input,
                    &mut sequence_store,
                    kmer_size,
                )?;
            info!(
                "Graph has {} nodes and {} edges",
                genome_graph.node_count(),
                genome_graph.edge_count()
            );

            if subcommand.linear_reduction {
                perform_linear_reduction(&mut genome_graph)?;
            }
            let genome_graph = genome_graph;

            info!("Computing maximal omnitigs");
            let mut omnitigs = Omnitigs::compute(&genome_graph);
            omnitigs.remove_reverse_complements(&genome_graph);
            print_omnitigs_statistics(&omnitigs, &mut latex_file)?;
            let mut omnitigs: Vec<_> = omnitigs.into_iter().map(Into::into).collect();
            if subcommand.linear_reduction {
                split_walks_at_node(
                    &mut omnitigs,
                    genome_graph.node_indices().last().unwrap(),
                    &genome_graph,
                )?;
                remove_subwalks_and_reverse_complements_from_walks(&mut omnitigs, &genome_graph);
            }

            info!(
                "Storing maximal omnitigs as fasta to '{}'",
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

            let mut genome_graph: PetWtdbg2Graph =
                genome_graph::io::wtdbg2::read_graph_from_wtdbg2_from_files(
                    nodes_file, reads_file, dot_file,
                )?;
            info!(
                "Graph has {} nodes and {} edges",
                genome_graph.node_count(),
                genome_graph.edge_count()
            );

            if subcommand.linear_reduction {
                perform_linear_reduction(&mut genome_graph)?;
            }
            let genome_graph = genome_graph;

            if !is_strongly_connected(&genome_graph) {
                warn!("Genome graph is not strongly connected!");
                warn!(
                    "Found {} SCCs",
                    decompose_strongly_connected_components(&genome_graph).len()
                );
            }

            info!("Computing maximal omnitigs");
            let mut omnitigs = Omnitigs::compute(&genome_graph);
            omnitigs.remove_reverse_complements(&genome_graph);
            print_omnitigs_statistics(&omnitigs, &mut latex_file)?;
            let mut omnitigs: Vec<_> = omnitigs.into_iter().map(Into::into).collect();
            if subcommand.linear_reduction {
                split_walks_at_node(
                    &mut omnitigs,
                    genome_graph.node_indices().last().unwrap(),
                    &genome_graph,
                )?;
                remove_subwalks_and_reverse_complements_from_walks(&mut omnitigs, &genome_graph);
            }

            if subcommand.output_as_wtdbg2_node_ids {
                info!("Storing omnitigs as node ids to '{}'", subcommand.output);
                genome_graph::io::wtdbg2::write_contigs_as_wtdbg2_node_ids_to_file(
                    &genome_graph,
                    omnitigs.iter(),
                    &subcommand.output,
                )?;
            } else {
                info!("Storing omnitigs as .ctg.lay to '{}'", subcommand.output);
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

            let mut genome_graph: PetWtdbg2DotGraph =
                genome_graph::io::wtdbg2::dot::read_graph_from_wtdbg2_dot_from_file(dot_file)?;

            info!(
                "Graph has {} nodes and {} edges",
                genome_graph.node_count(),
                genome_graph.edge_count()
            );

            if subcommand.linear_reduction {
                perform_linear_reduction(&mut genome_graph)?;
            }
            let genome_graph = genome_graph;

            info!("Computing maximal omnitigs");
            let mut omnitigs = Omnitigs::compute(&genome_graph);
            omnitigs.remove_reverse_complements(&genome_graph);
            print_omnitigs_statistics(&omnitigs, &mut latex_file)?;
            let mut omnitigs: Vec<_> = omnitigs.into_iter().map(Into::into).collect();
            if subcommand.linear_reduction {
                split_walks_at_node(
                    &mut omnitigs,
                    genome_graph.node_indices().last().unwrap(),
                    &genome_graph,
                )?;
                remove_subwalks_and_reverse_complements_from_walks(&mut omnitigs, &genome_graph);
            }

            info!("Storing omnitigs as node ids to '{}'", subcommand.output);
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
