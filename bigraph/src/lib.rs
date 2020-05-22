#[macro_use]
extern crate opaque_typedef_macros;

mod interface;
mod implementation;
mod index;

pub use interface::*;
pub use implementation::*;
pub use index::*;

#[cfg(test)]
mod tests {
    use crate::StaticGraph;

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

        fn print_graph<G: StaticGraph<i32, i32, usize> + std::fmt::Debug>(graph: &G) {
            println!("{:?}", graph);
        }

        print_graph(&graph);
    }
}
