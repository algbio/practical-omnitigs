use crate::CliOptions;
use clap::Clap;
use colored::*;
use genome_graph::bigraph::traitgraph::algo::components::{
    decompose_strongly_connected_components, decompose_weakly_connected_components,
    extract_subgraphs_from_node_mapping,
};
use genome_graph::bigraph::traitgraph::interface::{DynamicGraph, ImmutableGraphContainer};
use genome_graph::types::{PetBCalm2EdgeGraph, PetBCalm2NodeGraph};
use std::io::Write;
use omnitigs::macrotigs::macronodes::strongly_connected_macronode_algorithm::StronglyConnectedMacronodes;
use omnitigs::macrotigs::microtigs::strongly_connected_hydrostructure_based_maximal_microtig_algorithm::StronglyConnectedHydrostructureBasedMaximalMicrotigs;
use omnitigs::macrotigs::macronodes::MacronodeAlgorithm;
use omnitigs::macrotigs::microtigs::MaximalMicrotigsAlgorithm;
use omnitigs::macrotigs::macrotigs::default_macrotig_link_algorithm::DefaultMacrotigLinkAlgorithm;
use omnitigs::macrotigs::macrotigs::MaximalMacrotigsAlgorithm;
use traitsequence::interface::Sequence;

#[derive(Clap)]
pub struct VerifyEdgeCentricCommand {
    #[clap(
        short,
        long,
        about = "The kmer size selected when generating the input with bcalm2"
    )]
    pub kmer_size: usize,

    #[clap(
        short,
        long,
        about = "The output file, to which the graph should be written in bcalm2 format for verification purposes"
    )]
    pub output: Option<String>,

    #[clap(
        short,
        long,
        about = "A file to output the properties and statistics computed by this command formatted as a LaTeX table"
    )]
    pub latex: Option<String>,
}

#[derive(Clap)]
pub struct VerifyNodeCentricCommand {
    #[clap(
        short,
        long,
        about = "The output file, to which the graph should be written in bcalm2 format for verification purposes"
    )]
    pub output: Option<String>,
}

fn verify_components<Graph: Default + DynamicGraph, OutputWriter: std::io::Write>(
    genome_graph: &Graph,
    latex_file: &mut Option<OutputWriter>,
) -> crate::Result<()>
where
    Graph::NodeData: Clone,
    Graph::EdgeData: Clone,
{
    let wccs = decompose_weakly_connected_components(genome_graph);
    if wccs.is_empty() {
        warn!("No weakly connected components found");
    }

    let mut scc_amount_per_wcc = Vec::new();
    let mut scc_amount_per_wcc_float = Vec::new();
    let mut scc_count = 0;
    let mut non_scc_wcc_count = 0;
    for wcc in &wccs {
        let sccs =
            extract_subgraphs_from_node_mapping(wcc, &decompose_strongly_connected_components(wcc));
        scc_amount_per_wcc.push(sccs.len());
        scc_amount_per_wcc_float.push(sccs.len() as f64);
        scc_count += sccs.len();
        if sccs.len() >= 2 {
            non_scc_wcc_count += 1;
        }
    }
    let non_scc_wcc_string = if non_scc_wcc_count == 0 {
        "which all are strongly connected".to_string().normal()
    } else {
        format!("of which {} are not strongly connected", non_scc_wcc_count).red()
    };

    info!(
        "{} weakly connected components, {}",
        wccs.len(),
        non_scc_wcc_string
    );
    if let Some(latex_file) = latex_file.as_mut() {
        writeln!(latex_file, "WCCs & {} \\\\", wccs.len())?;
        writeln!(latex_file, "non-SCC WCCs & {} \\\\", non_scc_wcc_count)?;
    }

    let mut scc_amount_per_wcc_string = String::new();
    for &scc_amount in &scc_amount_per_wcc {
        if !scc_amount_per_wcc_string.is_empty() {
            scc_amount_per_wcc_string += ", ";
        }
        if scc_amount == 1 {
            scc_amount_per_wcc_string += "1";
        } else {
            assert!(scc_amount > 1);
            scc_amount_per_wcc_string += &format!("{}", scc_amount.to_string().yellow());
        }
    }

    info!("{} strongly connected components, that are distributed in the weakly connected components as [{}]", scc_count, scc_amount_per_wcc_string);
    if let Some(latex_file) = latex_file.as_mut() {
        writeln!(latex_file, "SCCs & {} \\\\", scc_count)?;
        if scc_amount_per_wcc.is_empty() {
            writeln!(latex_file, "max SCCs per WCC & N/A \\\\")?;
            writeln!(latex_file, "min SCCs per WCC & N/A \\\\")?;
            writeln!(latex_file, "median SCCs per WCC & N/A \\\\")?;
            writeln!(latex_file, "mean SCCs per WCC & N/A \\\\")?;
        } else {
            writeln!(
                latex_file,
                "max SCCs per WCC & {} \\\\",
                scc_amount_per_wcc.iter().max().unwrap()
            )?;
            writeln!(
                latex_file,
                "min SCCs per WCC & {} \\\\",
                scc_amount_per_wcc.iter().min().unwrap()
            )?;
            writeln!(
                latex_file,
                "median SCCs per WCC & {} \\\\",
                statistical::median(&scc_amount_per_wcc_float)
            )?;
            writeln!(
                latex_file,
                "mean SCCs per WCC & {} \\\\",
                statistical::mean(&scc_amount_per_wcc_float)
            )?;
        }
    }

    Ok(())
}

fn examine_macrotigs<Graph: Default + DynamicGraph, OutputWriter: std::io::Write>(
    genome_graph: &Graph,
    latex_file: &mut Option<OutputWriter>,
) -> crate::Result<()> {
    info!("");
    info!("=== Macrotigs ===");
    info!("");

    let macronodes = StronglyConnectedMacronodes::compute_macronodes(genome_graph);
    info!("{} macronodes", macronodes.len());
    if let Some(latex_file) = latex_file {
        writeln!(latex_file, "\\# macronodes & {} \\\\", macronodes.len())?;
    }

    let maximal_microtigs =
        StronglyConnectedHydrostructureBasedMaximalMicrotigs::compute_maximal_microtigs(
            genome_graph,
            &macronodes,
        );
    let macronodes_with_microtig_amount =
        macronodes.len() - maximal_microtigs.macronodes_without_microtig_amount();
    let macronodes_with_microtog_percent =
        100.0 * macronodes_with_microtig_amount as f64 / macronodes.len() as f64;
    info!(
        "{:.1}% of macronodes have at least one microtig",
        macronodes_with_microtog_percent
    );
    if let Some(latex_file) = latex_file {
        writeln!(
            latex_file,
            "\\% macronodes with at least one microtig & {:.1} \\\\",
            macronodes_with_microtog_percent
        )?;
    }

    let maximal_macrotigs =
        DefaultMacrotigLinkAlgorithm::compute_maximal_macrotigs(genome_graph, &maximal_microtigs);
    let total_macrotig_len: usize = maximal_macrotigs.iter().map(|m| m.len()).sum();
    info!("{} macrotigs", maximal_macrotigs.len());
    info!("total macrotig length: {}", total_macrotig_len);
    if let Some(latex_file) = latex_file {
        writeln!(
            latex_file,
            "\\# macrotigs & {} \\\\",
            maximal_macrotigs.len()
        )?;
        writeln!(
            latex_file,
            "total macrotig length & {} \\\\",
            total_macrotig_len
        )?;
    }

    if !maximal_macrotigs.is_empty() {
        let max_macrotig_len = maximal_macrotigs.iter().map(|m| m.len()).max().unwrap();
        let min_macrotig_len = maximal_macrotigs.iter().map(|m| m.len()).min().unwrap();
        let median_macrotig_len = statistical::median(
            &maximal_macrotigs
                .iter()
                .map(|m| m.len())
                .collect::<Vec<_>>(),
        );
        let mean_macrotig_len = statistical::mean(
            &maximal_macrotigs
                .iter()
                .map(|m| m.len() as f64)
                .collect::<Vec<_>>(),
        );

        info!("max macrotig edge length: {}", max_macrotig_len);
        info!("min macrotig edge length: {}", min_macrotig_len);
        info!("median macrotig edge length: {}", median_macrotig_len);
        info!("mean macrotig edge length: {:.1}", mean_macrotig_len);
        if let Some(latex_file) = latex_file {
            writeln!(
                latex_file,
                "max macrotig edge length & {} \\\\",
                max_macrotig_len
            )?;
            writeln!(
                latex_file,
                "min macrotig edge length & {} \\\\",
                min_macrotig_len
            )?;
            writeln!(
                latex_file,
                "median macrotig edge length & {} \\\\",
                median_macrotig_len
            )?;
            writeln!(
                latex_file,
                "mean macrotig edge length & {:.1} \\\\",
                mean_macrotig_len
            )?;
        }
    }

    Ok(())
}

pub(crate) fn verify_edge_centric(
    options: &CliOptions,
    subcommand: &VerifyEdgeCentricCommand,
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

    info!("");
    info!("========================");
    info!("=== Graph Statistics ===");
    info!("========================");
    info!("");
    info!(
        "{} nodes, {} edges",
        genome_graph.node_count(),
        genome_graph.edge_count()
    );
    if let Some(latex_file) = latex_file.as_mut() {
        writeln!(latex_file, "nodes & {} \\\\", genome_graph.node_count())?;
        writeln!(latex_file, "edges & {} \\\\", genome_graph.edge_count())?;
    }

    // Uncompacted unitigs
    let uncompacted_unitig_amount =
        omnitigs::unitigs::uncompacted_unitigs::count_uncompacted_edge_unitigs(&genome_graph);

    let log_string = format!("{} uncompacted unitigs", uncompacted_unitig_amount);
    if uncompacted_unitig_amount == 0 {
        info!("{}", log_string);
    } else {
        warn!("{}", log_string.yellow());
    }
    if let Some(latex_file) = latex_file.as_mut() {
        writeln!(
            latex_file,
            "uncompacted unitigs & {} \\\\",
            uncompacted_unitig_amount
        )?;
    }

    // Components
    verify_components(&genome_graph, &mut latex_file)?;

    // Macrotigs
    examine_macrotigs(&genome_graph, &mut latex_file)?;

    info!("");

    if let Some(output) = &subcommand.output {
        info!("Writing the unmodified bigraph to {}", output);
        genome_graph::io::bcalm2::write_edge_centric_bigraph_to_bcalm2_to_file(
            &genome_graph,
            output,
        )?;
    }
    Ok(())
}

pub(crate) fn verify_node_centric(
    options: &CliOptions,
    subcommand: &VerifyNodeCentricCommand,
) -> crate::Result<()> {
    info!("Reading bigraph from {}", options.input);
    let genome_graph: PetBCalm2NodeGraph =
        genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_node_centric_from_file(
            &options.input,
        )?;

    info!("");
    info!("========================");
    info!("=== Graph Statistics ===");
    info!("========================");
    info!("");
    info!(
        "{} binodes, {} edges",
        genome_graph.node_count() / 2,
        genome_graph.edge_count()
    );
    if genome_graph.node_count() % 2 != 0 {
        return Err(format!(
            "Internal error: uneven amount of nodes in a bidirected graph: {}",
            genome_graph.node_count()
        )
        .into());
    }

    // Uncompacted unitigs
    let uncompacted_unitig_amount =
        omnitigs::unitigs::uncompacted_unitigs::count_uncompacted_node_unitigs(&genome_graph);

    let log_string = format!("{} uncompacted unitigs", uncompacted_unitig_amount.total());
    let detail_string = format!(
        ", of which {} have length 2, {} have length 3 and {} are longer than 3",
        uncompacted_unitig_amount.len_2,
        uncompacted_unitig_amount.len_3,
        uncompacted_unitig_amount.len_4_more
    );
    if uncompacted_unitig_amount.total() == 0 {
        info!("{}", log_string);
    } else {
        warn!("{}{}", log_string.yellow(), detail_string.yellow());
    }

    /*if uncompacted_unitig_amount.total() % 2 != 0 {
        return Err(format!(
            "Internal error: uneven amount of total uncompacted unitigs in a bidirected graph: {}",
            uncompacted_unitig_amount.total()
        )
            .into());
    }
    if uncompacted_unitig_amount.len_2 % 2 != 0 {
        return Err(format!(
            "Internal error: uneven amount of len 2 uncompacted unitigs in a bidirected graph: {}",
            uncompacted_unitig_amount.len_2
        )
            .into());
    }
    if uncompacted_unitig_amount.len_3 % 2 != 0 {
        return Err(format!(
            "Internal error: uneven amount of len 3 uncompacted unitigs in a bidirected graph: {}",
            uncompacted_unitig_amount.len_3
        )
            .into());
    }
    if uncompacted_unitig_amount.len_4_more % 2 != 0 {
        return Err(format!("Internal error: uneven amount of len longer than 3 uncompacted unitigs in a bidirected graph: {}", uncompacted_unitig_amount.len_4_more).into());
    }*/

    // Components
    let mut latex_file: Option<std::fs::File> = None;
    verify_components(&genome_graph, &mut latex_file)?;

    info!("");

    if let Some(output) = &subcommand.output {
        info!("Writing the unmodified bigraph to {}", output);
        genome_graph::io::bcalm2::write_node_centric_bigraph_to_bcalm2_to_file(
            &genome_graph,
            output,
        )?;
    }
    Ok(())
}
