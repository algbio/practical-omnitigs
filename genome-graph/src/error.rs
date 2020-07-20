error_chain! {
    links {
        BCalm2IOError(crate::io::bcalm2::Error, crate::io::bcalm2::ErrorKind);
    }
}
