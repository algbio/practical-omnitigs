use crate::bigraph::interface::dynamic_bigraph::DynamicEdgeCentricBigraph;
use crate::bigraph::interface::dynamic_bigraph::DynamicNodeCentricBigraph;
use bigraph::interface::{dynamic_bigraph::DynamicBigraph, BidirectedData};
use bigraph::traitgraph::index::GraphIndex;
use bigraph::traitgraph::interface::GraphBase;
use bio::io::fasta::Record;
use compact_genome::{implementation::vector_genome_impl::VectorGenome, interface::Genome};
use num_traits::NumCast;
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use std::fmt::{Debug, Write};
use std::hash::Hash;
use std::iter::FromIterator;
use std::path::Path;

error_chain! {
    foreign_links {
        Io(std::io::Error);
        Fmt(std::fmt::Error);
    }

    errors {
        BCalm2IDError(id: String) {
            description("invalid node id")
            display("invalid node id: '{:?}'", id)
        }

        BCalm2LengthError(length: usize, sequence_length: usize) {
            description("the length in the description of a node does not match the length of its sequence")
            display("the length in the description of a node ({}) does not match the length of its sequence {}", length, sequence_length)
        }

        BCalm2UnknownParameterError(parameter: String) {
            description("unknown parameter")
            display("unknown parameter: '{:?}'", parameter)
        }

        BCalm2DuplicateParameterError(parameter: String) {
            description("duplicate parameter")
            display("duplicate parameter: '{:?}'", parameter)
        }

        BCalm2MalformedParameterError(parameter: String) {
            description("malformed parameter")
            display("malformed parameter: '{:?}'", parameter)
        }

        BCalm2MissingParameterError(parameter: String) {
            description("missing parameter")
            display("missing parameter: '{:?}'", parameter)
        }

        BCalm2NodeIdOutOfPrintingRange {
            description("node id is out of range (usize) for displaying")
            display("node id is out of range (usize) for displaying")
        }

        BCalm2NodeIdOutOfRange {
            description("node id is out of range (usize)")
            display("node id is out of range (usize)")
        }

        BCalm2NodeWithoutMirror {
            description("node has no mirror")
            display("node has no mirror")
        }

        BCalm2EdgeWithoutMirror {
            description("edge has no mirror")
            display("edge has no mirror")
        }
    }
}

#[derive(Debug)]
pub struct BCalm2NodeData {
    // TODO
}

/// The raw node data of a bcalm2 node, including edge information and redundant information (sequence length).
#[derive(Debug, Clone)]
pub struct PlainBCalm2NodeData {
    id: usize,
    sequence: VectorGenome,
    length: usize,
    total_abundance: usize,
    mean_abundance: f64,
    edges: Vec<PlainBCalm2Edge>,
}

/// The raw edge information of a bcalm2 node.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PlainBCalm2Edge {
    /// `true` means `+`, `false` means `-´
    from_side: bool,
    to_node: usize,
    /// `true` means `+`, `false` means `-´
    to_side: bool,
}

impl BidirectedData for PlainBCalm2NodeData {
    fn reverse_complement(&self) -> Self {
        let mut result = self.clone();
        result.sequence = result.sequence.reverse_complement();
        result
    }
}

/*impl PartialEq for PlainBCalm2NodeData {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmd::Ordering> {
        self.sequence.partial_cmp(other.sequence)
    }
}*/

impl PartialEq for PlainBCalm2NodeData {
    fn eq(&self, other: &Self) -> bool {
        self.sequence.eq(&other.sequence)
    }
}

impl Eq for PlainBCalm2NodeData {}

impl TryFrom<bio::io::fasta::Record> for PlainBCalm2NodeData {
    type Error = crate::error::Error;

    fn try_from(value: Record) -> crate::error::Result<Self> {
        let id = value
            .id()
            .parse()
            .map_err(|e| Error::with_chain(e, ErrorKind::BCalm2IDError(value.id().to_owned())))?;
        let sequence = VectorGenome::from_iter(value.seq()); // TODO store with bio
                                                             // TODO check if genome is valid

        let mut length = None;
        let mut total_abundance = None;
        let mut mean_abundance = None;
        let mut edges = Vec::new();

        for parameter in value.desc().unwrap_or("").split_whitespace() {
            ensure!(
                parameter.len() >= 5,
                Error::from(ErrorKind::BCalm2UnknownParameterError(parameter.to_owned()))
            );
            match &parameter[0..5] {
                "LN:i:" => {
                    ensure!(
                        length.is_none(),
                        Error::from(ErrorKind::BCalm2DuplicateParameterError(
                            parameter.to_owned()
                        ))
                    );
                    length = Some(parameter[5..].parse().map_err(|e| {
                        Error::with_chain(
                            e,
                            ErrorKind::BCalm2MalformedParameterError(parameter.to_owned()),
                        )
                    }));
                }
                "KC:i:" => {
                    ensure!(
                        total_abundance.is_none(),
                        Error::from(ErrorKind::BCalm2DuplicateParameterError(
                            parameter.to_owned()
                        ))
                    );
                    total_abundance = Some(parameter[5..].parse().map_err(|e| {
                        Error::with_chain(
                            e,
                            ErrorKind::BCalm2MalformedParameterError(parameter.to_owned()),
                        )
                    }));
                }
                "KM:f:" | "km:f:" => {
                    ensure!(
                        mean_abundance.is_none(),
                        Error::from(ErrorKind::BCalm2DuplicateParameterError(
                            parameter.to_owned()
                        ))
                    );
                    mean_abundance = Some(parameter[5..].parse().map_err(|e| {
                        Error::with_chain(
                            e,
                            ErrorKind::BCalm2MalformedParameterError(parameter.to_owned()),
                        )
                    }));
                }
                _ => match &parameter[0..2] {
                    "L:" => {
                        let parts: Vec<_> = parameter.split(':').collect();
                        ensure!(
                            parts.len() == 4,
                            Error::from(ErrorKind::BCalm2MalformedParameterError(
                                parameter.to_owned()
                            ))
                        );
                        let forward_reverse_to_bool = |c| match c {
                            "+" => Ok(true),
                            "-" => Ok(false),
                            _ => Err(Error::from(ErrorKind::BCalm2MalformedParameterError(
                                parameter.to_owned(),
                            ))),
                        };
                        let from_side = forward_reverse_to_bool(parts[1])?;
                        let to_node = parts[2].parse().map_err(|e| {
                            Error::with_chain(
                                e,
                                ErrorKind::BCalm2MalformedParameterError(parameter.to_owned()),
                            )
                        })?;
                        let to_side = forward_reverse_to_bool(parts[3])?;
                        edges.push(PlainBCalm2Edge {
                            from_side,
                            to_node,
                            to_side,
                        });
                    }
                    _ => bail!(Error::from(ErrorKind::BCalm2UnknownParameterError(
                        parameter.to_owned()
                    ))),
                },
            }
        }

        let length = length.unwrap_or_else(|| {
            bail!(Error::from(ErrorKind::BCalm2MissingParameterError(
                "length (LN)".to_owned()
            )))
        })?;
        ensure!(
            length == sequence.len(),
            Error::from(ErrorKind::BCalm2LengthError(length, sequence.len()))
        );
        let total_abundance = total_abundance.unwrap_or_else(|| {
            bail!(Error::from(ErrorKind::BCalm2MissingParameterError(
                "total abundance (KC)".to_owned()
            )))
        })?;
        let mean_abundance = mean_abundance.unwrap_or_else(|| {
            bail!(Error::from(ErrorKind::BCalm2MissingParameterError(
                "mean abundance (KM)".to_owned()
            )))
        })?;

        Ok(Self {
            id,
            sequence,
            length,
            total_abundance,
            mean_abundance,
            edges,
        })
    }
}

impl<'a> From<&'a PlainBCalm2NodeData> for PlainBCalm2NodeData {
    fn from(data: &'a PlainBCalm2NodeData) -> Self {
        data.clone()
    }
}

/////////////////////////////
////// NODE CENTRIC IO //////
/////////////////////////////

pub fn read_bigraph_from_bcalm2_as_node_centric_from_file<
    P: AsRef<Path>,
    NodeData: From<PlainBCalm2NodeData> + BidirectedData,
    EdgeData: Default + Clone,
    Graph: DynamicNodeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    path: P,
) -> crate::error::Result<Graph> {
    read_bigraph_from_bcalm2_as_node_centric(
        bio::io::fasta::Reader::from_file(path).map_err(Error::from)?,
    )
}

pub fn read_bigraph_from_bcalm2_as_node_centric<
    R: std::io::Read,
    NodeData: From<PlainBCalm2NodeData> + BidirectedData,
    EdgeData: Default + Clone,
    Graph: DynamicNodeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    reader: bio::io::fasta::Reader<R>,
) -> crate::error::Result<Graph> {
    struct BiEdge {
        from_node: usize,
        plain_edge: PlainBCalm2Edge,
    }

    let mut bigraph = Graph::default();
    let mut edges = Vec::new();

    for record in reader.records() {
        let record: PlainBCalm2NodeData = record.map_err(Error::from)?.try_into()?;
        edges.extend(record.edges.iter().map(|e| BiEdge {
            from_node: record.id,
            plain_edge: e.clone(),
        }));
        let record_id = record.id;
        let id = bigraph.add_node(record.into());
        assert_eq!(id, record_id.into());
    }

    bigraph.add_mirror_nodes();
    assert!(bigraph.verify_node_pairing());

    for edge in edges {
        let from_node = if edge.plain_edge.from_side {
            edge.from_node.into()
        } else {
            bigraph.mirror_node(edge.from_node.into()).unwrap()
        };
        let to_node = if edge.plain_edge.to_side {
            edge.plain_edge.to_node.into()
        } else {
            bigraph.mirror_node(edge.plain_edge.to_node.into()).unwrap()
        };
        bigraph.add_edge(from_node, to_node, EdgeData::default());
    }

    bigraph.add_node_centric_mirror_edges();
    assert!(bigraph.verify_node_mirror_property());
    Ok(bigraph)
}

fn write_plain_bcalm2_node_data_to_bcalm2(
    node: &PlainBCalm2NodeData,
    out_neighbors: Vec<(bool, usize, bool)>,
) -> crate::error::Result<String> {
    let mut result = String::new();
    write!(
        result,
        "LN:i:{} KC:i:{} km:f:{:.1}",
        node.length, node.total_abundance, node.mean_abundance
    )
    .map_err(Error::from)?;
    for (node_type, neighbor_id, neighbor_type) in out_neighbors {
        write!(
            result,
            " L:{}:{}:{}",
            if node_type { "+" } else { "-" },
            <usize as NumCast>::from(neighbor_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfPrintingRange))?,
            if neighbor_type { "+" } else { "-" }
        )
        .map_err(Error::from)?;
    }
    Ok(result)
}

pub fn write_node_centric_bigraph_to_bcalm2_to_file<
    P: AsRef<Path>,
    NodeData, //: Into<PlainBCalm2NodeData<IndexType>>,
    EdgeData: Default + Clone,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    graph: &Graph,
    path: P,
) -> crate::error::Result<()>
where
    for<'a> PlainBCalm2NodeData: From<&'a NodeData>,
{
    write_node_centric_bigraph_to_bcalm2(
        graph,
        bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

pub fn write_node_centric_bigraph_to_bcalm2<
    W: std::io::Write,
    NodeData,
    EdgeData: Default + Clone,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    graph: &Graph,
    mut writer: bio::io::fasta::Writer<W>,
) -> crate::error::Result<()>
where
    for<'a> PlainBCalm2NodeData: From<&'a NodeData>,
{
    let mut output_nodes = vec![false; graph.node_count()];

    for node_id in graph.node_indices() {
        if !output_nodes[graph
            .mirror_node(node_id)
            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutMirror))?
            .as_usize()]
        {
            output_nodes[node_id.as_usize()] = true;
        }
    }

    for node_id in graph.node_indices() {
        if output_nodes[node_id.as_usize()] {
            let node_data = PlainBCalm2NodeData::from(graph.node_data(node_id));
            let mirror_node_id = graph
                .mirror_node(node_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutMirror))?;
            /*let mirror_node_data = PlainBCalm2NodeData::<IndexType>::from(
                graph
                    .node_data(mirror_node_id)
                    .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?,
            );*/
            let mut out_neighbors_plus = Vec::new();
            let mut out_neighbors_minus = Vec::new();

            for neighbor in graph.out_neighbors(node_id) {
                let neighbor_node_id = neighbor.node_id.as_usize();

                out_neighbors_plus.push((
                    true,
                    if output_nodes[neighbor_node_id] {
                        neighbor.node_id.as_usize()
                    } else {
                        graph
                            .mirror_node(neighbor.node_id)
                            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?
                            .as_usize()
                    },
                    output_nodes[neighbor_node_id],
                ));
            }
            for neighbor in graph.out_neighbors(mirror_node_id) {
                let neighbor_node_id = neighbor.node_id.as_usize();

                out_neighbors_minus.push((
                    false,
                    if output_nodes[neighbor_node_id] {
                        neighbor.node_id.as_usize()
                    } else {
                        graph
                            .mirror_node(neighbor.node_id)
                            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?
                            .as_usize()
                    },
                    output_nodes[neighbor_node_id],
                ));
            }

            out_neighbors_plus.sort();
            out_neighbors_minus.sort();
            out_neighbors_plus.append(&mut out_neighbors_minus);
            let out_neighbors = out_neighbors_plus;

            let mut printed_node_id = String::new();
            write!(printed_node_id, "{}", node_data.id).map_err(Error::from)?;
            let node_description =
                write_plain_bcalm2_node_data_to_bcalm2(&node_data, out_neighbors)?;
            let node_sequence = node_data.sequence.into_vec();

            writer
                .write(&printed_node_id, Some(&node_description), &node_sequence)
                .map_err(Error::from)?;
        }
    }

    Ok(())
}

/////////////////////////////
////// EDGE CENTRIC IO //////
/////////////////////////////

pub fn read_bigraph_from_bcalm2_as_edge_centric_from_file<
    P: AsRef<Path>,
    NodeData: Default + Clone,
    EdgeData: From<PlainBCalm2NodeData> + Clone + Eq + BidirectedData,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    path: P,
    kmer_size: usize,
) -> crate::error::Result<Graph> {
    read_bigraph_from_bcalm2_as_edge_centric(
        bio::io::fasta::Reader::from_file(path).map_err(Error::from)?,
        kmer_size,
    )
}

fn get_or_create_node<Graph: DynamicBigraph, G: Genome + Hash>(
    bigraph: &mut Graph,
    id_map: &mut HashMap<G, <Graph as GraphBase>::NodeIndex>,
    genome: &G,
) -> <Graph as GraphBase>::NodeIndex
where
    for<'a> &'a G: IntoIterator<Item = u8>,
    <Graph as GraphBase>::NodeData: Default,
    <Graph as GraphBase>::EdgeData: Clone,
{
    if let Some(node) = id_map.get(genome) {
        *node
    } else {
        let node = bigraph.add_node(Default::default());
        id_map.insert(genome.clone(), node);

        let reverse_complement = genome.reverse_complement();
        if &reverse_complement == genome {
            bigraph.set_mirror_nodes(node, node);
        } else {
            let mirror_node = bigraph.add_node(Default::default());
            id_map.insert(reverse_complement, mirror_node);
            bigraph.set_mirror_nodes(node, mirror_node);
        }

        node
    }
}

pub fn read_bigraph_from_bcalm2_as_edge_centric<
    R: std::io::Read,
    NodeData: Default + Clone,
    EdgeData: From<PlainBCalm2NodeData> + Clone + Eq + BidirectedData,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    reader: bio::io::fasta::Reader<R>,
    kmer_size: usize,
) -> crate::error::Result<Graph>
where
    <Graph as GraphBase>::NodeIndex: Clone,
{
    let mut bigraph = Graph::default();
    let mut id_map = HashMap::new();
    let node_kmer_size = kmer_size - 1;

    for record in reader.records() {
        let record: PlainBCalm2NodeData = record.map_err(Error::from)?.try_into()?;
        let reverse_complement = record.reverse_complement();

        let pre_plus = record.sequence.prefix(node_kmer_size);
        let pre_minus = reverse_complement.sequence.prefix(node_kmer_size);
        let succ_plus = record.sequence.suffix(node_kmer_size);
        let succ_minus = reverse_complement.sequence.suffix(node_kmer_size);

        let pre_plus = get_or_create_node(&mut bigraph, &mut id_map, &pre_plus);
        let pre_minus = get_or_create_node(&mut bigraph, &mut id_map, &pre_minus);
        let succ_plus = get_or_create_node(&mut bigraph, &mut id_map, &succ_plus);
        let succ_minus = get_or_create_node(&mut bigraph, &mut id_map, &succ_minus);

        bigraph.add_edge(pre_plus, succ_plus, record.clone().into());
        bigraph.add_edge(pre_minus, succ_minus, record.reverse_complement().into());
    }

    assert!(bigraph.verify_node_pairing());
    assert!(bigraph.verify_edge_mirror_property());
    Ok(bigraph)
}

pub fn write_edge_centric_bigraph_to_bcalm2_to_file<
    P: AsRef<Path>,
    NodeData, //: Into<PlainBCalm2NodeData<IndexType>>,
    EdgeData: BidirectedData + Clone + Eq,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    graph: &Graph,
    path: P,
) -> crate::error::Result<()>
where
    for<'a> PlainBCalm2NodeData: From<&'a EdgeData>,
{
    write_edge_centric_bigraph_to_bcalm2(
        graph,
        bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

pub fn write_edge_centric_bigraph_to_bcalm2<
    W: std::io::Write,
    NodeData,
    EdgeData: BidirectedData + Clone + Eq,
    Graph: DynamicEdgeCentricBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    graph: &Graph,
    mut writer: bio::io::fasta::Writer<W>,
) -> crate::error::Result<()>
where
    for<'a> PlainBCalm2NodeData: From<&'a EdgeData>,
{
    let mut output_edges = vec![false; graph.edge_count()];

    for edge_id in graph.edge_indices() {
        if !output_edges[graph
            .mirror_edge_edge_centric(edge_id)
            .ok_or_else(|| Error::from(ErrorKind::BCalm2EdgeWithoutMirror))?
            .as_usize()]
        {
            output_edges[edge_id.as_usize()] = true;
        }
    }

    for edge_id in graph.edge_indices() {
        if output_edges[edge_id.as_usize()] {
            let node_data = PlainBCalm2NodeData::from(graph.edge_data(edge_id));
            let mirror_edge_id = graph
                .mirror_edge_edge_centric(edge_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2EdgeWithoutMirror))?;
            let to_node_plus = graph.edge_endpoints(edge_id).to_node;
            let to_node_minus = graph.edge_endpoints(mirror_edge_id).to_node;

            let mut out_neighbors_plus = Vec::new();
            let mut out_neighbors_minus = Vec::new();

            for neighbor in graph.out_neighbors(to_node_plus) {
                let neighbor_edge_id = neighbor.edge_id.as_usize();

                out_neighbors_plus.push((
                    true,
                    if output_edges[neighbor_edge_id] {
                        PlainBCalm2NodeData::from(graph.edge_data(neighbor.edge_id)).id
                    } else {
                        PlainBCalm2NodeData::from(
                            graph.edge_data(
                                graph
                                    .mirror_edge_edge_centric(neighbor.edge_id)
                                    .ok_or_else(|| {
                                        Error::from(ErrorKind::BCalm2EdgeWithoutMirror)
                                    })?,
                            ),
                        )
                        .id
                    },
                    output_edges[neighbor_edge_id],
                ));
            }
            for neighbor in graph.out_neighbors(to_node_minus) {
                let neighbor_edge_id = neighbor.edge_id.as_usize();

                out_neighbors_minus.push((
                    false,
                    if output_edges[neighbor_edge_id] {
                        PlainBCalm2NodeData::from(graph.edge_data(neighbor.edge_id)).id
                    } else {
                        PlainBCalm2NodeData::from(
                            graph.edge_data(
                                graph
                                    .mirror_edge_edge_centric(neighbor.edge_id)
                                    .ok_or_else(|| {
                                        Error::from(ErrorKind::BCalm2EdgeWithoutMirror)
                                    })?,
                            ),
                        )
                        .id
                    },
                    output_edges[neighbor_edge_id],
                ));
            }

            out_neighbors_plus.sort();
            out_neighbors_minus.sort();
            out_neighbors_plus.append(&mut out_neighbors_minus);
            let out_neighbors = out_neighbors_plus;

            let mut printed_node_id = String::new();
            write!(printed_node_id, "{}", node_data.id).map_err(Error::from)?;
            let node_description =
                write_plain_bcalm2_node_data_to_bcalm2(&node_data, out_neighbors)?;
            let node_sequence = node_data.sequence.into_vec();

            writer
                .write(&printed_node_id, Some(&node_description), &node_sequence)
                .map_err(Error::from)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::io::bcalm2::{
        read_bigraph_from_bcalm2_as_edge_centric, read_bigraph_from_bcalm2_as_node_centric,
        write_edge_centric_bigraph_to_bcalm2, write_node_centric_bigraph_to_bcalm2,
    };
    use crate::types::{PetBCalm2EdgeGraph, PetBCalm2NodeGraph};

    #[test]
    fn test_node_read_write() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:1:-\n\
            AGT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            GGTCTCGGGTAAGT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:1:-\n\
            ATGATG\n";
        let input = Vec::from(test_file);

        let graph: PetBCalm2NodeGraph =
            read_bigraph_from_bcalm2_as_node_centric(bio::io::fasta::Reader::new(test_file))
                .unwrap();
        let mut output = Vec::new();
        write_node_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:1:-\n\
            AGT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            AATCTCGGGTAAAC\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:1:-\n\
            ACGAGG\n";
        let input = Vec::from(test_file);

        let graph: PetBCalm2EdgeGraph =
            read_bigraph_from_bcalm2_as_edge_centric(bio::io::fasta::Reader::new(test_file), 3)
                .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_self_loops() {
        let test_file: &'static [u8] =
            b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:0:+ L:+:1:- L:-:0:- L:-:2:+\n\
            AAA\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            GGTCTCGGGTAATT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:0:+ L:-:1:-\n\
            TTGATG\n";
        let input = Vec::from(test_file);
        println!("{}", String::from_utf8(input.clone()).unwrap());

        let graph: PetBCalm2EdgeGraph =
            read_bigraph_from_bcalm2_as_edge_centric(bio::io::fasta::Reader::new(test_file), 3)
                .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_plus_minus_loop() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:0:- L:+:1:- L:+:2:+\n\
            CAT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:1:- L:+:2:+\n\
            GGTCTCGGGTAAAT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:0:- L:-:1:- L:-:2:+\n\
            ATGATT\n";
        let input = Vec::from(test_file);
        println!("{}", String::from_utf8(input.clone()).unwrap());

        let graph: PetBCalm2EdgeGraph =
            read_bigraph_from_bcalm2_as_edge_centric(bio::io::fasta::Reader::new(test_file), 3)
                .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_minus_plus_loop() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:1:- L:-:0:+\n\
            ATG\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            GGTCTCGGGTAACA\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:1:-\n\
            CAGATT\n";
        let input = Vec::from(test_file);
        println!("{}", String::from_utf8(input.clone()).unwrap());

        let graph: PetBCalm2EdgeGraph =
            read_bigraph_from_bcalm2_as_edge_centric(bio::io::fasta::Reader::new(test_file), 3)
                .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_plus_minus_and_minus_plus_loop() {
        let test_file: &'static [u8] =
            b">0 LN:i:4 KC:i:4 km:f:3.0 L:+:0:- L:+:1:- L:+:2:+ L:-:0:+\n\
            CGAT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:1:- L:+:2:+\n\
            GGTCTCGGGTAAAT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:0:- L:-:1:- L:-:2:+\n\
            ATGATG\n";
        let input = Vec::from(test_file);
        println!("{}", String::from_utf8(input.clone()).unwrap());

        let graph: PetBCalm2EdgeGraph =
            read_bigraph_from_bcalm2_as_edge_centric(bio::io::fasta::Reader::new(test_file), 3)
                .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_all_loops() {
        let test_file: &'static [u8] = b">0 LN:i:4 KC:i:4 km:f:3.0 L:+:0:- L:+:0:+ L:+:1:- L:+:2:+ L:-:0:- L:-:0:+ L:-:1:- L:-:2:+\n\
            ATAT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:0:+ L:+:1:- L:+:2:+\n\
            TCTCGGAAGTAAAT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:0:- L:-:0:+ L:-:1:- L:-:2:+\n\
            ATGATG\n";
        let input = Vec::from(test_file);

        let graph: PetBCalm2EdgeGraph =
            read_bigraph_from_bcalm2_as_edge_centric(bio::io::fasta::Reader::new(test_file), 3)
                .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_forward_merge() {
        let test_file: &'static [u8] = b"\
            >0 LN:i:3 KC:i:4 km:f:3.0 L:+:2:-\n\
            AGT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:2:-\n\
            GGTCTCGGGTAAGT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:+:0:- L:+:1:-\n\
            AAGAAC\n";
        let input = Vec::from(test_file);

        let graph: PetBCalm2EdgeGraph =
            read_bigraph_from_bcalm2_as_edge_centric(bio::io::fasta::Reader::new(test_file), 3)
                .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }

    #[test]
    fn test_edge_read_write_multigraph() {
        let test_file: &'static [u8] = b"\
            >0 LN:i:7 KC:i:4 km:f:3.0 L:+:2:+ L:-:2:-\n\
            AGTTCTC\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:2:+ L:-:2:-\n\
            AGTCTCGGGTAATC\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:+:0:+ L:+:1:+ L:-:0:- L:-:1:-\n\
            TCGAAG\n";
        let input = Vec::from(test_file);

        let graph: PetBCalm2EdgeGraph =
            read_bigraph_from_bcalm2_as_edge_centric(bio::io::fasta::Reader::new(test_file), 3)
                .unwrap();
        let mut output = Vec::new();
        write_edge_centric_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output))
            .unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }
}
