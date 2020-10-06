use crate::index::GraphIndex;
use crate::interface::{DynamicGraph, StaticGraph};
use std::io::{BufRead, BufReader, Read, Write};

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

/// Read a graph as Hamiltonian circuit problem encoded as TSP in TSPLIB format (used by concorde).
pub fn read_hamcircuit_from_concorde_tsp<Graph: DynamicGraph, Reader: Read>(
    graph: &mut Graph,
    reader: &mut Reader,
) where
    Graph::NodeData: Default,
    Graph::EdgeData: Default,
{
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

    let mut node_count = None;

    let mut reader = BufReader::new(reader);
    loop {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();

        if line.trim().starts_with("DIMENSION:") {
            let line = line.trim().split(' ').last().unwrap();
            node_count = Some(line.parse::<usize>().unwrap());
        }

        if line.trim().starts_with("EDGE_WEIGHT_SECTION") {
            break;
        }
    }

    let node_count = node_count.unwrap();

    for _ in 0..node_count {
        graph.add_node(Default::default());
    }

    for n1 in graph.node_indices() {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        let edges = line.trim().split(' ');

        for (n2, edge_weight) in graph.node_indices().zip(edges) {
            let edge_weight: usize = edge_weight.trim().parse().unwrap();
            if edge_weight == 1 {
                graph.add_edge(n1, n2, Default::default());
            }
        }
    }

    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    assert_eq!(line.trim(), "EOF");
}
