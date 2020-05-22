pub trait Genome<'a>: From<&'a [u8]> + Into<Vec<u8>> + std::fmt::Display {
    fn reverse_complement(&self) -> Self;
}

fn ascii_complement(char: u8) -> Option<u8> {
    match char {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'G' => Some(b'C'),
        b'C' => Some(b'G'),
        _ => None,
    }
}