use crate::interface::static_bigraph::StaticEdgeCentricBigraph;
use crate::interface::BidirectedData;
use crate::traitgraph::index::GraphIndex;
use bitvector::BitVector;
use std::fmt::Write;

/// Algorithms related to bidirected Eulerian graphs.
pub mod eulerian;
/// Algorithms related to covering a bigraph with biwalks.
pub mod walk_cover;

fn bitvector_to_index_string(bitvector: &BitVector) -> String {
    let mut result = String::new();
    for i in 0..bitvector.capacity() {
        if bitvector.contains(i) {
            write!(result, "{} ", i).unwrap()
        }
    }
    result
}

fn mark_edge_and_mirror<
    EdgeData: BidirectedData + Eq,
    Graph: StaticEdgeCentricBigraph<EdgeData = EdgeData>,
>(
    bitvector: &mut BitVector,
    graph: &Graph,
    edge_index: Graph::EdgeIndex,
) {
    debug_assert!(!bitvector.contains(edge_index.as_usize()));
    bitvector.insert(edge_index.as_usize());
    //println!("Marked edge {}", edge_index.as_usize());

    let mirror_edge = graph.mirror_edge_edge_centric(edge_index).unwrap();
    if bitvector.insert(mirror_edge.as_usize()) {
        //println!("Marked edge {}", mirror_edge.as_usize());
    } else {
        panic!(
            "bitvector {}\nforwards {:?}\nmirrors {:?}",
            bitvector_to_index_string(bitvector),
            graph
                .edges_between(
                    graph.edge_endpoints(edge_index).from_node,
                    graph.edge_endpoints(edge_index).to_node
                )
                .map(|e| e.as_usize())
                .collect::<Vec<_>>(),
            graph.topological_mirror_edges(edge_index)
        );
    }
}
