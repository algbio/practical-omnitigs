use crate::CliOptions;
use clap::Clap;
use colored::*;
use genome_graph::bigraph::traitgraph::algo::components::{
    decompose_strongly_connected_components, decompose_weakly_connected_components,
    extract_subgraphs_from_node_mapping,
};
use genome_graph::bigraph::traitgraph::interface::{DynamicGraph, ImmutableGraphContainer};
use genome_graph::types::{PetBCalm2EdgeGraph, PetBCalm2NodeGraph};

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

fn verify_components<Graph: Default + DynamicGraph>(genome_graph: &Graph)
where
    Graph::NodeData: Clone,
    Graph::EdgeData: Clone,
{
    let wccs = decompose_weakly_connected_components(genome_graph);
    let mut scc_amount_per_wcc = Vec::new();
    let mut scc_count = 0;
    let mut non_scc_wcc_count = 0;
    for wcc in &wccs {
        let sccs =
            extract_subgraphs_from_node_mapping(wcc, &decompose_strongly_connected_components(wcc));
        scc_amount_per_wcc.push(sccs.len());
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

    let mut scc_amount_per_wcc_string = String::new();
    for scc_amount in scc_amount_per_wcc {
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
}

pub(crate) fn verify_edge_centric(
    options: &CliOptions,
    subcommand: &VerifyEdgeCentricCommand,
) -> crate::Result<()> {
    let genome_graph: PetBCalm2EdgeGraph =
        genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file(
            &options.input,
            subcommand.kmer_size,
        )?;

    info!(
        "Reading bigraph from '{}' with kmer size {}",
        options.input, subcommand.kmer_size
    );
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

    // Uncompacted unitigs
    let uncompacted_unitig_amount =
        omnitigs::unitigs::uncompacted_unitigs::count_uncompacted_edge_unitigs(&genome_graph);

    let log_string = format!("{} uncompacted unitigs", uncompacted_unitig_amount);
    if uncompacted_unitig_amount == 0 {
        info!("{}", log_string);
    } else {
        warn!("{}", log_string.yellow());
    }

    // Components
    verify_components(&genome_graph);

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
    verify_components(&genome_graph);

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
