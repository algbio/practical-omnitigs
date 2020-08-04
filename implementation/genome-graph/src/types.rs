pub type PetBCalm2Graph = crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
    crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
        crate::io::bcalm2::PlainBCalm2NodeData,
        (),
        usize,
    >,
>;
