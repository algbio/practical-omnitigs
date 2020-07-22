pub type PetBCalm2Graph = crate::bigraph::NodeBigraphWrapper<
    crate::io::bcalm2::PlainBCalm2NodeData<usize>,
    (),
    usize,
    crate::bigraph::petgraph::Graph<
        crate::io::bcalm2::PlainBCalm2NodeData<usize>,
        (),
        crate::bigraph::petgraph::Directed,
        usize,
    >,
>;
