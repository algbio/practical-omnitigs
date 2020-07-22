#![recursion_limit = "1024"]
#[macro_use]
extern crate error_chain;

mod error;
pub mod io;
pub mod types;

pub use bigraph;
pub use error::*;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
