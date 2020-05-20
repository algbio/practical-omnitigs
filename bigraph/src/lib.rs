#[macro_use]
extern crate opaque_typedef_macros;

mod bigraph;

pub use bigraph::{EdgeIndex, NodeIndex, ImmutableGraphContainer, MutableGraphContainer, Bigraph};

#[cfg(test)]
mod tests {
    use crate::bigraph::StaticGraph;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn test_type_layout() {
        let mut graph = petgraph::Graph::default();
        let n1 = graph.add_node(4);
        let n2 = graph.add_node(5);
        graph.add_edge(n1, n2, 6);

        let graph_box = Box::new(graph);
        let _static_graph_box: Box<dyn StaticGraph<i32, i32, usize>> = graph_box;
    }
}
