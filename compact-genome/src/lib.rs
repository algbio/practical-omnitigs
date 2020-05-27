mod genome;
mod vector_genome_impl;

pub use genome::*;
pub use vector_genome_impl::*;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
