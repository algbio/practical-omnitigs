#[macro_use]
extern crate opaque_typedef_macros;

mod bigraph;

pub use bigraph::{EdgeIndex, NodeIndex, StaticGraph, DynamicGraph, Bigraph};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
