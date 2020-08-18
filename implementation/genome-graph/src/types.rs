pub type PetBCalm2NodeGraph =
    crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            crate::io::bcalm2::PlainBCalm2NodeData,
            (),
            usize,
        >,
    >;

pub type PetBCalm2EdgeGraph =
    crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
        crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
            (),
            crate::io::bcalm2::PlainBCalm2NodeData,
            usize,
        >,
    >;
