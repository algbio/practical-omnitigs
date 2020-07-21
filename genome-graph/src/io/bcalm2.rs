use bigraph::{BidirectedNodeData, DynamicBigraph};
use bio::io::fasta::Record;
use num_traits::PrimInt;
use std::convert::{TryFrom, TryInto};
use std::fmt::Debug;
use std::path::Path;
use std::str::FromStr;

error_chain! {
    foreign_links {
        Io(std::io::Error);
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
    }
}

#[derive(Debug)]
pub struct BCalm2NodeData {
    // TODO
}

/// The raw node data of a bcalm2 node, including edge information and redundant information (sequence length).
#[derive(Debug, Clone)]
pub struct PlainBCalm2NodeData<IndexType> {
    id: IndexType,
    sequence: Vec<u8>,
    length: usize,
    total_abundance: usize,
    mean_abundance: f64,
    edges: Vec<PlainBCalm2Edge<IndexType>>,
}

/// The raw edge information of a bcalm2 node.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PlainBCalm2Edge<IndexType> {
    /// `true` means `+`, `false` means `-´
    from_side: bool,
    to_node: IndexType,
    /// `true` means `+`, `false` means `-´
    to_side: bool,
}

impl<IndexType> BidirectedNodeData for PlainBCalm2NodeData<IndexType> {
    fn reverse_complement(&self) -> Self {
        // let mut result = self.clone();
        // result.sequence = // TODO bio reverse complement
        unimplemented!()
    }
}

impl<IndexType: FromStr> TryFrom<bio::io::fasta::Record> for PlainBCalm2NodeData<IndexType>
where
    <IndexType as FromStr>::Err: std::error::Error + Send + 'static,
{
    type Error = crate::Error;

    fn try_from(value: Record) -> crate::Result<Self> {
        let id = value
            .id()
            .parse()
            .map_err(|e| Error::with_chain(e, ErrorKind::BCalm2IDError(value.id().to_owned())))?;
        let sequence = value.seq().to_owned(); // TODO store with bio

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

pub fn load_bigraph_from_bcalm2<
    P: AsRef<Path>,
    NodeData: From<PlainBCalm2NodeData<IndexType>>,
    EdgeData: Default,
    IndexType: PrimInt + FromStr + Debug,
    T: DynamicBigraph<NodeData, EdgeData, IndexType> + Default,
>(
    path: P,
) -> crate::Result<T>
where
    <IndexType as FromStr>::Err: std::error::Error + Send + 'static,
{
    struct BiEdge<IndexType> {
        _from_node: IndexType,
        _plain_edge: PlainBCalm2Edge<IndexType>,
    }

    let mut digraph = T::default();
    let mut edges = Vec::new();

    for record in bio::io::fasta::Reader::from_file(path)
        .map_err(Error::from)?
        .records()
    {
        let record: PlainBCalm2NodeData<IndexType> = record.map_err(Error::from)?.try_into()?;
        edges.extend(record.edges.iter().map(|e| BiEdge {
            _from_node: record.id,
            _plain_edge: e.clone(),
        }));
        let record_id = record.id;
        let id = digraph.add_node(record.into());
        assert_eq!(id, record_id.into());
    }

    for _edge in edges {
        unimplemented!()
        // TODO digraph.add_edge()
    }

    unimplemented!()
}
