use std::path::Path;
use std::fs::File;
use std::io::{BufReader, Lines, BufRead};
use crate::Error;

pub fn open_input_text_file<P: AsRef<Path>>(path: P) -> Result<Lines<BufReader<File>>, Error> {
    let file = File::open(path)?;
    Ok(BufReader::new(file).lines())
}