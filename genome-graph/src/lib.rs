mod error;
mod io;

pub use bigraph;
pub use error::*;
pub use io::*;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
