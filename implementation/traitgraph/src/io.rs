use crate::index::GraphIndex;
use crate::interface::StaticGraph;
use std::io::Write;

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

/// Write the graph as Hamiltonian circuit problem encoded as TSP in TSPLIB format (used by concorde).
pub fn write_hamcircuit_as_concorde_tsp<Graph: StaticGraph, Writer: Write>(
    graph: &Graph,
    writer: &mut Writer,
) {
    // NAME:  demo
    // TYPE: TSP
    // COMMENT: 4 vertexes asymmetric problem
    // DIMENSION:  4
    // EDGE_WEIGHT_TYPE: EXPLICIT
    // EDGE_WEIGHT_FORMAT: FULL_MATRIX
    // EDGE_WEIGHT_SECTION
    // 9999 11 8 4
    // 10 9999 7 2
    // 6 5 9999 4
    // 6 3 9 9999
    // EOF

    writeln!(writer, "NAME: none\nTYPE: TSP\nCOMMENT: none\nDIMENSION: {}\nEDGE_WEIGHT_TYPE: EXPLICIT\nEDGE_WEIGHT_FORMAT: FULL_MATRIX\nEDGE_WEIGHT_SECTION", graph.node_count()).unwrap();

    for n1 in graph.node_indices() {
        for n2 in graph.node_indices() {
            if n1 == n2 || !graph.contains_edge_between(n1, n2) {
                write!(writer, " {}", graph.node_count() * 10).unwrap();
            } else {
                write!(writer, " 1").unwrap();
            }
        }

        writeln!(writer).unwrap();
    }

    writeln!(writer, "EOF").unwrap();
}
