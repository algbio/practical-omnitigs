error_chain! {
    foreign_links {
        // For some weird reasons I don't understand, the doc comments have to be put after the item in this macro...
        Io(std::io::Error)
        /// An IO error.
        ;
        Anyhow(anyhow::Error)
        /// Any error passed through anyhow.
        ;
    }

    links {
        // For some weird reasons I don't understand, the doc comments have to be put after the item in this macro...
        BCalm2IoError(crate::io::bcalm2::Error, crate::io::bcalm2::ErrorKind)
        /// A wrapper for errors thrown by bcalm2 IO.
        ;
        FastaIoError(crate::io::fasta::Error, crate::io::fasta::ErrorKind)
        /// A wrapper for errors thrown by fasta IO.
        ;
    }

    errors {
        #[allow(missing_docs)]
        GfaUnknownOverlapPattern {
            description("an L-line was encountered, but the overlap pattern is unknown")
            display("an L-line was encountered, but the overlap pattern is unknown")
        }

        #[allow(missing_docs)]
        GfaMissingOverlapPattern {
            description("an L-line was encountered, but the overlap pattern is missing")
            display("an L-line was encountered, but the overlap pattern is missing")
        }

        #[allow(missing_docs)]
        GfaMissingNode {
            description("an L-line was encountered, at least one of the nodes is missing")
            display("an L-line was encountered, at least one of the nodes is missing")
        }
    }
}
