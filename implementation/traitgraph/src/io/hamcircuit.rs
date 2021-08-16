use crate::interface::{DynamicGraph, StaticGraph};
use std::io::{BufRead, BufReader, Read, Write};

/// Write the graph as Hamiltonian circuit problem encoded as ATSP in TSPLIB format (used by concorde).
pub fn write_hamcircuit_as_tsplib_atsp<Graph: StaticGraph, Writer: Write>(
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

    writeln!(writer, "NAME: none\nTYPE: ATSP\nCOMMENT: none\nDIMENSION: {}\nEDGE_WEIGHT_TYPE: EXPLICIT\nEDGE_WEIGHT_FORMAT: FULL_MATRIX\nEDGE_WEIGHT_SECTION", graph.node_count()).unwrap();

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

/// Read a graph as Hamiltonian circuit problem encoded as ATSP in TSPLIB format (used by concorde).
pub fn read_hamcircuit_from_tsplib_atsp<Graph: DynamicGraph, Reader: Read>(
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
    debug_assert_eq!(line.trim(), "EOF");
}

/// Write the graph as Hamiltonian circuit problem encoded as TSP in TSPLIB format (used by concorde).
pub fn write_hamcircuit_as_tsplib_tsp<Graph: StaticGraph, Writer: Write>(
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

    writeln!(writer, "NAME: none\nTYPE: TSP\nCOMMENT: none\nDIMENSION: {}\nEDGE_WEIGHT_TYPE: EXPLICIT\nEDGE_WEIGHT_FORMAT: FULL_MATRIX\nEDGE_WEIGHT_SECTION", graph.node_count() * 2).unwrap();

    // First half of the matrix.
    for n2 in graph.node_indices() {
        for _ in graph.node_indices() {
            write!(writer, " {}", graph.node_count() * 10).unwrap();
        }

        for n1 in graph.node_indices() {
            if n1 == n2 {
                write!(writer, " 0").unwrap();
            } else if graph.contains_edge_between(n1, n2) {
                write!(writer, " 10").unwrap();
            } else {
                write!(writer, " 12").unwrap();
            }
        }

        writeln!(writer).unwrap();
    }

    // Second half of the matrix.
    for n1 in graph.node_indices() {
        for n2 in graph.node_indices() {
            if n1 == n2 {
                write!(writer, " 0").unwrap();
            } else if graph.contains_edge_between(n1, n2) {
                write!(writer, " 10").unwrap();
            } else {
                write!(writer, " 12").unwrap();
            }
        }

        for _ in graph.node_indices() {
            write!(writer, " {}", graph.node_count() * 10).unwrap();
        }

        writeln!(writer).unwrap();
    }

    writeln!(writer, "EOF").unwrap();
}

/// Read a graph as Hamiltonian circuit problem encoded as TSP in TSPLIB format (used by concorde).
pub fn read_hamcircuit_from_tsplib_tsp<Graph: DynamicGraph, Reader: Read>(
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
            debug_assert_eq!(node_count.unwrap() % 2, 0, "Node count is odd");
            node_count = Some(node_count.unwrap() / 2);
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

        for (n2, edge_weight) in graph.node_indices().zip(edges.skip(node_count)) {
            let edge_weight: usize = edge_weight.trim().parse().unwrap();
            if edge_weight == 10 {
                graph.add_edge(n2, n1, Default::default());
            }
        }
    }

    for _ in graph.node_indices() {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
    }

    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    debug_assert_eq!(line.trim(), "EOF");
}

#[cfg(test)]
mod tests {
    use crate::implementation::petgraph_impl;
    use crate::interface::{ImmutableGraphContainer, MutableGraphContainer};
    use crate::io::hamcircuit::{read_hamcircuit_from_tsplib_tsp, write_hamcircuit_as_tsplib_tsp};
    use std::io::{BufReader, BufWriter};

    #[test]
    fn test_write_read_simple() {
        let mut graph = petgraph_impl::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n0, n2, ());
        graph.add_edge(n1, n2, ());
        graph.add_edge(n1, n3, ());
        graph.add_edge(n2, n1, ());
        graph.add_edge(n3, n0, ());

        let mut writer = BufWriter::new(Vec::new());
        write_hamcircuit_as_tsplib_tsp(&graph, &mut writer);
        let mut result = petgraph_impl::new::<(), ()>();
        let buffer = writer.into_inner().unwrap();
        let mut reader = BufReader::new(buffer.as_slice());
        read_hamcircuit_from_tsplib_tsp(&mut result, &mut reader);
        debug_assert_eq!(graph.node_count(), result.node_count());
        debug_assert_eq!(graph.edge_count(), result.edge_count());

        for n1 in graph.node_indices() {
            for n2 in graph.node_indices() {
                debug_assert_eq!(
                    graph.contains_edge_between(n1, n2),
                    result.contains_edge_between(n1, n2)
                );
            }
        }
    }
}
