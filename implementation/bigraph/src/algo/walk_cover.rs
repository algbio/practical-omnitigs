use crate::algo::mark_edge_and_mirror;
use crate::interface::static_bigraph::StaticEdgeCentricBigraph;
use crate::interface::BidirectedData;
use crate::traitgraph::index::GraphIndex;
use bitvector::BitVector;
use traitgraph::walks::VecEdgeWalk;

/// Computes a set of biwalks covering the graph, without any further guarantees.
pub fn arbitrary_biwalk_cover<
    EdgeData: BidirectedData + Eq,
    Graph: StaticEdgeCentricBigraph<EdgeData = EdgeData>,
>(
    graph: &Graph,
) -> Vec<VecEdgeWalk<Graph>> {
    let mut used_edges = BitVector::new(graph.edge_count());
    let mut walks = Vec::new();

    for edge_index in graph.edge_indices() {
        if used_edges.contains(edge_index.as_usize()) {
            continue;
        }
        //println!("Starting new cycle with edge {}", edge_index.as_usize());

        let mut walk = Vec::new();

        //println!("Start edge {}", start_edge_index.as_usize());
        // Extend backwards first
        mark_edge_and_mirror(&mut used_edges, graph, edge_index);
        walk.push(edge_index);
        let mut current_node = graph.edge_endpoints(edge_index).from_node;

        let mut has_neighbor = true;
        while has_neighbor {
            //println!("Expanding node {}", current_node.as_usize());
            has_neighbor = false;
            for neighbor in graph.in_neighbors(current_node) {
                if !used_edges.contains(neighbor.edge_id.as_usize()) {
                    walk.push(neighbor.edge_id);
                    mark_edge_and_mirror(&mut used_edges, graph, neighbor.edge_id);
                    current_node = neighbor.node_id;
                    has_neighbor = true;
                    break;
                }
            }
        }

        walk.reverse();

        // Extend forwards
        has_neighbor = true;
        current_node = graph.edge_endpoints(edge_index).to_node;
        while has_neighbor {
            //println!("Expanding node {}", current_node.as_usize());
            has_neighbor = false;
            for neighbor in graph.out_neighbors(current_node) {
                if !used_edges.contains(neighbor.edge_id.as_usize()) {
                    walk.push(neighbor.edge_id);
                    mark_edge_and_mirror(&mut used_edges, graph, neighbor.edge_id);
                    current_node = neighbor.node_id;
                    has_neighbor = true;
                    break;
                }
            }
        }

        walks.push(VecEdgeWalk::new(walk));
    }

    walks
}
