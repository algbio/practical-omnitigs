#[derive(Debug)]
pub enum Error {
    IoError(std::io::Error),
}

impl From<std::io::Error> for Error {
    fn from(io_error: std::io::Error) -> Self {
        Self::IoError(io_error)
    }
}
