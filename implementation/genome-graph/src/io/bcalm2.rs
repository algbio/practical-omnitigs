use bigraph::interface::{dynamic_bigraph::DynamicBigraph, BidirectedNodeData};
use bigraph::traitgraph::index::GraphIndex;
use bio::io::fasta::Record;
use compact_genome::{implementation::vector_genome_impl::VectorGenome, interface::Genome};
use num_traits::NumCast;
use std::convert::{TryFrom, TryInto};
use std::fmt::{Debug, Write};
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

        BCalm2NodeWithoutPartner {
            description("node has no partner")
            display("node has no partner")
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

impl BidirectedNodeData for PlainBCalm2NodeData {
    fn reverse_complement(&self) -> Self {
        let mut result = self.clone();
        result.sequence = result.sequence.reverse_complement();
        result
    }
}

impl TryFrom<bio::io::fasta::Record> for PlainBCalm2NodeData {
    type Error = crate::error::Error;

    fn try_from(value: Record) -> crate::error::Result<Self> {
        let id = value
            .id()
            .parse()
            .map_err(|e| Error::with_chain(e, ErrorKind::BCalm2IDError(value.id().to_owned())))?;
        let sequence = VectorGenome::from_iter(value.seq()); // TODO store with bio

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

pub fn read_bigraph_from_bcalm2_file<
    P: AsRef<Path>,
    NodeData: From<PlainBCalm2NodeData>,
    EdgeData: Default + Clone,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    path: P,
) -> crate::error::Result<Graph> {
    read_bigraph_from_bcalm2(bio::io::fasta::Reader::from_file(path).map_err(Error::from)?)
}

pub fn read_bigraph_from_bcalm2<
    R: std::io::Read,
    NodeData: From<PlainBCalm2NodeData>,
    EdgeData: Default + Clone,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
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

    bigraph.add_partner_nodes();
    assert!(bigraph.verify_node_pairing());

    for edge in edges {
        let from_node = if edge.plain_edge.from_side {
            edge.from_node.into()
        } else {
            bigraph.partner_node(edge.from_node.into()).unwrap()
        };
        let to_node = if edge.plain_edge.to_side {
            edge.plain_edge.to_node.into()
        } else {
            bigraph
                .partner_node(edge.plain_edge.to_node.into())
                .unwrap()
        };
        bigraph.add_edge(from_node, to_node, EdgeData::default());
    }

    bigraph.add_mirror_edges();
    assert!(bigraph.verify_mirror_property());
    Ok(bigraph)
}

fn write_binode_to_bcalm2(
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

pub fn write_bigraph_to_bcalm2_file<
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
    write_bigraph_to_bcalm2(
        graph,
        bio::io::fasta::Writer::to_file(path).map_err(Error::from)?,
    )
}

pub fn write_bigraph_to_bcalm2<
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
            .partner_node(node_id)
            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutPartner))?
            .as_usize()]
        {
            output_nodes[node_id.as_usize()] = true;
        }
    }

    for node_id in graph.node_indices() {
        if output_nodes[node_id.as_usize()] {
            let node_data = PlainBCalm2NodeData::from(graph.node_data(node_id));
            let partner_node_id = graph
                .partner_node(node_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutPartner))?;
            /*let partner_node_data = PlainBCalm2NodeData::<IndexType>::from(
                graph
                    .node_data(partner_node_id)
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
                            .partner_node(neighbor.node_id)
                            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?
                            .as_usize()
                    },
                    output_nodes[neighbor_node_id],
                ));
            }
            for neighbor in graph.out_neighbors(partner_node_id) {
                let neighbor_node_id = neighbor.node_id.as_usize();

                out_neighbors_minus.push((
                    false,
                    if output_nodes[neighbor_node_id] {
                        neighbor.node_id.as_usize()
                    } else {
                        graph
                            .partner_node(neighbor.node_id)
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
            let node_description = write_binode_to_bcalm2(&node_data, out_neighbors)?;
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
    use crate::io::bcalm2::{read_bigraph_from_bcalm2, write_bigraph_to_bcalm2};
    use crate::types::PetBCalm2Graph;

    #[test]
    fn test_read_write() {
        let test_file: &'static [u8] = b">0 LN:i:3 KC:i:4 km:f:3.0 L:+:1:-\n\
            AGT\n\
            >1 LN:i:14 KC:i:2 km:f:3.2 L:+:0:- L:+:2:+\n\
            GGTCTCGGGTAAGT\n\
            >2 LN:i:6 KC:i:15 km:f:2.2 L:-:1:-\n\
            ATGATG\n";
        let input = Vec::from(test_file);

        let graph: PetBCalm2Graph =
            read_bigraph_from_bcalm2(bio::io::fasta::Reader::new(test_file)).unwrap();
        let mut output = Vec::new();
        write_bigraph_to_bcalm2(&graph, bio::io::fasta::Writer::new(&mut output)).unwrap();

        assert_eq!(
            input,
            output,
            "in:\n{}\n\nout:\n{}\n",
            String::from_utf8(input.clone()).unwrap(),
            String::from_utf8(output.clone()).unwrap()
        );
    }
}
