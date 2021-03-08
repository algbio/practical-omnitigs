use crate::error::Result;
use bigraph::interface::dynamic_bigraph::DynamicBigraph;
use compact_genome::implementation::vector_genome_impl::VectorGenome;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::iter::FromIterator;
use std::path::Path;
use std::rc::Rc;

/// Type of graphs read from gfa files.
pub type PetGFAGraph<NodeData, EdgeData> = crate::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper<
    crate::bigraph::traitgraph::implementation::petgraph_impl::petgraph::graph::DiGraph<
        BidirectedGFANodeData<NodeData>,
        EdgeData,
        usize,
    >,
>;

/// Node data of a bidirected graph read from GFA
pub struct BidirectedGFANodeData<T> {
    /// The sequence of this node. If forward is false, then this must be reverse complemented.
    pub sequence: Rc<VectorGenome>,
    /// True if this node is the forward node of sequence, false if it is the reverse complement node.
    pub forward: bool,
    /// Further data.
    pub data: T,
}

/// Read an edge-centric bigraph in gfa format from a file.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_edge_centric_bigraph_from_file<
    P: AsRef<Path>,
    NodeData: Default,
    EdgeData: Default,
    Graph: DynamicBigraph<NodeData = BidirectedGFANodeData<NodeData>, EdgeData = EdgeData> + Default,
>(
    gfa_file: P,
) -> Result<(Graph, usize)> {
    read_gfa_as_edge_centric_bigraph(BufReader::new(File::open(gfa_file)?))
}

/// Read an edge-centric bigraph in gfa format from a `BufRead`.
/// This method also returns the k-mer length given in the gfa file.
pub fn read_gfa_as_edge_centric_bigraph<
    R: BufRead,
    NodeData: Default,
    EdgeData: Default,
    Graph: DynamicBigraph<NodeData = BidirectedGFANodeData<NodeData>, EdgeData = EdgeData> + Default,
>(
    gfa: R,
) -> Result<(Graph, usize)> {
    let mut graph = Graph::default();
    let mut k = usize::max_value();

    for line in gfa.lines() {
        let line = line?;

        if line.starts_with('H') {
            assert!(graph.is_empty());
            for column in line.split('\t') {
                if let Some(stripped) = column.strip_prefix("KL:Z:") {
                    assert_eq!(k, usize::max_value());
                    k = stripped.parse().unwrap();
                }
            }
        } else if line.starts_with('S') {
            assert_eq!(graph.edge_count(), 0);
            assert_ne!(k, usize::max_value());

            let mut columns = line.split('\t').skip(1);
            let node_index: usize = columns.next().unwrap().parse().unwrap();
            assert_eq!((node_index - 1) * 2, graph.node_count());

            let sequence = Rc::new(VectorGenome::from_iter(columns.next().unwrap().bytes()));
            let n1 = graph.add_node(BidirectedGFANodeData {
                sequence: sequence.clone(),
                forward: true,
                data: Default::default(),
            });
            let n2 = graph.add_node(BidirectedGFANodeData {
                sequence: sequence.clone(),
                forward: false,
                data: Default::default(),
            });
            graph.set_mirror_nodes(n1, n2);
        } else if line.starts_with('L') {
            assert_ne!(k, usize::max_value());

            let mut columns = line.split('\t').skip(1);
            let n1 = (columns.next().unwrap().parse::<usize>().unwrap()
                + if columns.next().unwrap() == "+" { 0 } else { 1 })
            .into();
            let n2 = (columns.next().unwrap().parse::<usize>().unwrap()
                + if columns.next().unwrap() == "+" { 0 } else { 1 })
            .into();

            graph.add_edge(n1, n2, Default::default());
            graph.add_edge(
                graph.mirror_node(n2).unwrap(),
                graph.mirror_node(n1).unwrap(),
                Default::default(),
            );
        }
    }

    Ok((graph, k))
}
