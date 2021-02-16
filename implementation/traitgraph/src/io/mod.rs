use crate::index::GraphIndex;
use crate::interface::StaticGraph;
use std::io::Write;

/// IO methods related to hamiltonian circuits and the concorde TSP solver.
pub mod hamcircuit;

/// Write the graph in the following format, ignoring node and edge data.
///
/// ```text
/// <node count> <edge count>
/// <from node> <to node>
/// ```
///
/// The second line is repeated for each edge.
pub fn write_topology<Graph: StaticGraph, Writer: Write>(graph: &Graph, writer: &mut Writer) {
    writeln!(writer, "{} {}", graph.node_count(), graph.edge_count()).unwrap();
    for node in graph.node_indices() {
        for out_neighbor in graph.out_neighbors(node) {
            writeln!(
                writer,
                "{} {}",
                node.as_usize(),
                out_neighbor.node_id.as_usize()
            )
            .unwrap();
        }
    }
}
