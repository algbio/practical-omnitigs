use crate::interface::StaticGraph;

/// Returns true if the graph contains a Eulerian cycle.
pub fn decomposes_into_eulerian_cycles<Graph: StaticGraph>(graph: &Graph) -> bool {
    for node_index in graph.node_indices() {
        if graph.in_degree(node_index) != graph.out_degree(node_index) {
            return false;
        }
    }

    true
}

/// Compute a vector of nodes that has indegree != outdegree.
pub fn find_non_eulerian_nodes<Graph: StaticGraph>(graph: &Graph) -> Vec<Graph::NodeIndex> {
    let mut result = Vec::new();
    for node_index in graph.node_indices() {
        if graph.in_degree(node_index) != graph.out_degree(node_index) {
            result.push(node_index);
        }
    }
    result
}

/// Compute a vector of tuples of nodes and outdegree - indegree that has indegree != outdegree.
pub fn find_non_eulerian_nodes_with_differences<Graph: StaticGraph>(
    graph: &Graph,
) -> Vec<(Graph::NodeIndex, isize)> {
    let mut node_indices_and_differences = Vec::new();
    for node_index in graph.node_indices() {
        let difference =
            graph.out_degree(node_index) as isize - graph.in_degree(node_index) as isize;
        if difference != 0 {
            node_indices_and_differences.push((node_index, difference));
        }
    }
    node_indices_and_differences
}
