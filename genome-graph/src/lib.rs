mod error;
mod io;

pub use error::*;
pub use io::*;
pub use bigraph;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
