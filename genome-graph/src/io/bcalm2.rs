use crate::error::Error;
use bigraph::StaticBigraph;
use num_traits::PrimInt;
use std::path::Path;
use crate::io::util::open_input_text_file;

pub struct BCalm2NodeData {
    // TODO
}

pub fn load_bigraph_from_bcalm2<P: AsRef<Path>, NodeData: From<BCalm2NodeData>, EdgeData: Default, IndexType: PrimInt, T: StaticBigraph<NodeData, EdgeData, IndexType>>(path: P) -> Result<T, Error> {
    let mut lines = open_input_text_file(path)?;

    println!("{:?}", lines.next());

    unimplemented!()
}