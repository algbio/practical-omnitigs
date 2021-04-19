use crate::CliOptions;
use clap::Clap;
use genome_graph::bigraph::traitgraph::algo::components::is_strongly_connected;
use genome_graph::bigraph::traitgraph::algo::predefined_graphs::{
    compute_m_from_n_and_c, create_random_graph,
};
use genome_graph::bigraph::traitgraph::implementation::petgraph_impl;
use genome_graph::bigraph::traitgraph::io::hamcircuit::{
    read_hamcircuit_from_tsplib_tsp, write_hamcircuit_as_tsplib_tsp,
};
use omnitigs::hamiltonian::preprocess_hamiltonian_circuit;
use omnitigs::macrotigs::macrotigs::Macrotigs;
use omnitigs::node_covering_node_visible_one_circular_safe::compute_maximal_node_covering_node_visible_one_circular_safe_walks;
use omnitigs::traitgraph::interface::ImmutableGraphContainer;
use omnitigs::traitgraph::interface::MutableGraphContainer;
use omnitigs::traitgraph::walks::NodeWalk;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::time::Instant;
use traitsequence::interface::Sequence;

#[derive(Clap)]
pub struct HamCircuitCommand {
    #[clap(short, long, about = "The input file in tsplib format")]
    pub input: String,

    #[clap(
        short,
        long,
        about = "If given, generate random strongly connected graphs instead of reading an input. Must have an argument of the form n<node count>+c<arc factor>."
    )]
    random: Option<String>,

    #[clap(long, about = "The file to write the raw graph in TSPLIB format.")]
    output_raw: String,

    #[clap(
        long,
        about = "The file to write the preprocessed graph in TSPLIB format."
    )]
    output_preprocessed: String,
}

#[allow(clippy::unnecessary_wraps)]
pub(crate) fn hamcircuit(
    _options: &CliOptions,
    subcommand: &HamCircuitCommand,
) -> crate::Result<()> {
    let mut graph = petgraph_impl::new::<(), ()>();

    if let Some(random) = &subcommand.random {
        let (node_count, c_value) = scan_fmt!(random, "n{d}+c{f}", usize, f64).expect("Could not parse argument of random. Make sure it fulfills the format in the help message.");
        let mut loop_count = 0;

        loop {
            graph.clear();
            create_random_graph(&mut graph, node_count, c_value, &mut rand::thread_rng());
            if is_strongly_connected(&graph) {
                break;
            }
            loop_count += 1;
            if loop_count == 10 {
                info!("Did not find a strongly connected graph after 10 attempts, using n = {} and c = {} (which translates to an edge/node ration of {:.2}%)",
                      node_count,
                      c_value,
                      compute_m_from_n_and_c(node_count, c_value) as f64 / node_count as f64 * 100.0
                );
            }
        }

        info!(
            "Generated strongly connected graph with {} nodes and {} arcs",
            graph.node_count(),
            graph.edge_count()
        );
    } else {
        info!("Reading graph from input file: '{}'", &subcommand.input);
        let input_file = File::open(&subcommand.input).expect("Input file does not exist");
        let mut input_reader = BufReader::new(input_file);
        read_hamcircuit_from_tsplib_tsp(&mut graph, &mut input_reader);
    };

    info!("Opening raw output file: '{}'", subcommand.output_raw);
    let output_raw_file = File::create(&subcommand.output_raw).unwrap();
    let mut output_raw_writer = BufWriter::new(output_raw_file);
    info!(
        "Opening preprocessed output file: '{}'",
        subcommand.output_preprocessed
    );
    let output_preprocessed_file = File::create(&subcommand.output_preprocessed).unwrap();
    let mut output_preprocessed_writer = BufWriter::new(output_preprocessed_file);

    info!("Writing raw graph");
    write_hamcircuit_as_tsplib_tsp(&graph, &mut output_raw_writer);
    output_raw_writer.flush().unwrap();
    drop(output_raw_writer);

    let safe_walks = compute_maximal_node_covering_node_visible_one_circular_safe_walks(
        &graph,
        &Macrotigs::compute(&graph),
    );
    info!("Found {} safe walks", safe_walks.len());
    for walk in &safe_walks {
        info!("{:?}", walk);
    }
    info!(
        "Found {} safe walks of length > 2",
        safe_walks.iter().filter(|walk| walk.len() > 2).count()
    );
    info!(
        "Found {} non-trivial safe walks",
        safe_walks
            .iter()
            .filter(|walk| walk.is_non_trivial(&graph))
            .count()
    );

    info!("Preprocessing for the Hamiltonian circuit problem");
    let start_time = Instant::now();
    let preprocessed = preprocess_hamiltonian_circuit(&graph);
    let duration = Instant::now() - start_time;
    info!("Preprocessing took {:.1} seconds", duration.as_secs_f64());

    if let Some(preprocessed) = preprocessed {
        info!(
            "Preprocessing removed {:.1}% of nodes and {:.1}% of arcs",
            (graph.node_count() - preprocessed.node_count()) as f64 / graph.node_count() as f64
                * 100.0,
            (graph.edge_count() - preprocessed.edge_count()) as f64 / graph.edge_count() as f64
                * 100.0,
        );

        info!("Writing preprocessed graph");
        write_hamcircuit_as_tsplib_tsp(&preprocessed, &mut output_preprocessed_writer);
    } else {
        info!("Preprocessing the graph revealed that it is not hamiltonian");
        info!("Writing dummy graph");
        graph.clear();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n1, n2, ());
        graph.add_edge(n2, n1, ());
        graph.add_edge(n1, n0, ());
        write_hamcircuit_as_tsplib_tsp(&graph, &mut output_preprocessed_writer);
    }

    Ok(())
}
