error_chain! {
    links {
        // For some weird reasons I don't understand, the doc comments have to be put after the item in this macro...
        BCalm2IOError(crate::io::bcalm2::Error, crate::io::bcalm2::ErrorKind)
        /// A wrapper for errors thrown by bcalm2 IO.
        ;
        FastaIOError(crate::io::fasta::Error, crate::io::fasta::ErrorKind)
        /// A wrapper for errors thrown by fasta IO.
        ;
    }
}
