use bigraph::{BidirectedNodeData, DynamicBigraph, NodeIndex};
use bio::io::fasta::Record;
use compact_genome::{Genome, VectorGenome};
use num_traits::{NumCast, PrimInt};
use std::convert::{TryFrom, TryInto};
use std::fmt::{Debug, Display, Write};
use std::iter::FromIterator;
use std::path::Path;
use std::str::FromStr;

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
pub struct PlainBCalm2NodeData<IndexType: PrimInt> {
    id: IndexType,
    sequence: VectorGenome,
    length: usize,
    total_abundance: usize,
    mean_abundance: f64,
    edges: Vec<PlainBCalm2Edge<IndexType>>,
}

/// The raw edge information of a bcalm2 node.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PlainBCalm2Edge<IndexType: PrimInt> {
    /// `true` means `+`, `false` means `-´
    from_side: bool,
    to_node: IndexType,
    /// `true` means `+`, `false` means `-´
    to_side: bool,
}

impl<IndexType: PrimInt> BidirectedNodeData for PlainBCalm2NodeData<IndexType> {
    fn reverse_complement(&self) -> Self {
        let mut result = self.clone();
        result.sequence = result.sequence.reverse_complement();
        result
    }
}

impl<IndexType: FromStr + PrimInt> TryFrom<bio::io::fasta::Record>
    for PlainBCalm2NodeData<IndexType>
where
    <IndexType as FromStr>::Err: std::error::Error + Send + 'static,
{
    type Error = crate::Error;

    fn try_from(value: Record) -> crate::Result<Self> {
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
                // A bug in bcalm2 causes it to output the lower-case pattern, instead of the documented upper-case one.
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

impl<'a, IndexType: PrimInt> From<&'a PlainBCalm2NodeData<IndexType>> for PlainBCalm2NodeData<IndexType> {
    fn from(data: &'a PlainBCalm2NodeData<IndexType>) -> Self {
        data.clone()
    }
}

pub fn read_bigraph_from_bcalm2<
    P: AsRef<Path>,
    NodeData: From<PlainBCalm2NodeData<IndexType>>,
    EdgeData: Default + Clone,
    IndexType: PrimInt + FromStr + Debug,
    T: DynamicBigraph<NodeData, EdgeData, IndexType> + Default,
>(
    path: P,
) -> crate::Result<T>
where
    <IndexType as FromStr>::Err: std::error::Error + Send + 'static,
{
    struct BiEdge<IndexType: PrimInt> {
        from_node: IndexType,
        plain_edge: PlainBCalm2Edge<IndexType>,
    }

    let mut bigraph = T::default();
    let mut edges = Vec::new();

    for record in bio::io::fasta::Reader::from_file(path)
        .map_err(Error::from)?
        .records()
    {
        let record: PlainBCalm2NodeData<IndexType> = record.map_err(Error::from)?.try_into()?;
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

fn write_binode_to_bcalm2<IndexType: PrimInt + Debug>(
    node: &PlainBCalm2NodeData<IndexType>,
    out_neighbors: Vec<(NodeIndex<IndexType>, bool)>,
) -> crate::Result<String> {
    let mut result = String::new();
    write!(
        result,
        "LN:i:{} KC:i:{} km:f:{}",
        node.length, node.total_abundance, node.mean_abundance
    )
    .map_err(Error::from)?;
    for (neighbor_id, neighbor_type) in out_neighbors {
        write!(
            result,
            " L:+:{}:{}",
            <usize as NumCast>::from(neighbor_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfPrintingRange))?,
            if neighbor_type { "+" } else { "-" }
        )
        .map_err(Error::from)?;
    }
    Ok(result)
}

pub fn write_bigraph_to_bcalm2<
    P: AsRef<Path>,
    NodeData, //: Into<PlainBCalm2NodeData<IndexType>>,
    EdgeData: Default + Clone,
    IndexType: PrimInt + Debug + Display,
    T: DynamicBigraph<NodeData, EdgeData, IndexType> + Default,
>(
    graph: &T,
    path: P,
) -> crate::Result<()>
where
    PlainBCalm2NodeData<IndexType>: for<'a> From<&'a NodeData>,
{
    let mut output_nodes = vec![false; graph.node_count()];
    let mut writer = bio::io::fasta::Writer::to_file(path).map_err(Error::from)?;

    for node_id in graph.node_indices() {
        if !output_nodes[<usize as NumCast>::from(
            graph
                .partner_node(node_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeWithoutPartner))?,
        )
        .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?]
        {
            output_nodes[<usize as NumCast>::from(node_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?] = true;
        }
    }

    for node_id in graph.node_indices() {
        if output_nodes[<usize as NumCast>::from(node_id)
            .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?]
        {
            let node_data = PlainBCalm2NodeData::<IndexType>::from(
                graph
                    .node_data(node_id)
                    .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?,
            );
            let mut out_neighbors = Vec::new();
            for neighbor in graph
                .out_neighbors(node_id)
                .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?
            {
                out_neighbors.push((
                    neighbor.node_id,
                    output_nodes[<usize as NumCast>::from(neighbor.node_id)
                        .ok_or_else(|| Error::from(ErrorKind::BCalm2NodeIdOutOfRange))?],
                ));
            }
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
