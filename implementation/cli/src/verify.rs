use crate::CliOptions;
use colored::*;
use genome_graph::bigraph::traitgraph::interface::ImmutableGraphContainer;
use genome_graph::types::PetBCalm2Graph;

pub(crate) fn verify(options: &CliOptions) -> crate::Result<()> {
    info!("Reading bigraph from {}", options.input);
    let genome_graph: PetBCalm2Graph =
        genome_graph::io::bcalm2::read_bigraph_from_bcalm2_file(&options.input)?;

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
        genome_graph::bigraph::traitgraph::algo::unitigs::count_uncompacted_node_unitigs(
            &genome_graph,
        );

    let log_string = format!(
        "{} uncompacted unitigs",
        uncompacted_unitig_amount.total()
    );
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
    let wccs =
        genome_graph::bigraph::traitgraph::algo::components::decompose_weakly_connected_components(
            &genome_graph,
        );
    let mut non_scc_wcc_count = 0;
    for wcc in &wccs {
        if !genome_graph::bigraph::traitgraph::algo::components::is_strongly_connected(wcc) {
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

    info!("");

    if let Some(output) = &options.output {
        info!("Writing the unmodified bigraph to {}", output);
        genome_graph::io::bcalm2::write_bigraph_to_bcalm2_file(&genome_graph, output)?;
    }
    Ok(())
}
