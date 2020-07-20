use crate::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Lines};
use std::path::Path;

pub fn open_input_text_file<P: AsRef<Path>>(path: P) -> Result<Lines<BufReader<File>>, Error> {
    let file = File::open(path)?;
    Ok(BufReader::new(file).lines())
}
