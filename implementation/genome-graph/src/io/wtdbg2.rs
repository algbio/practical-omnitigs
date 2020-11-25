use crate::error::Result;
use bigraph::traitgraph::interface::DynamicGraph;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Node data as given in a .1.nodes file from wtdbg2.
pub struct PlainWtdbg2NodeData {
    /// The index of the node in wtdbg2.
    pub index: usize,
    /// The read associations of the node.
    pub read_associations: Vec<Wtdbg2ReadAssociation>,
}

/// Read associations of nodes of wtdbg2.
pub struct Wtdbg2ReadAssociation {
    /// The identifier of the read.
    pub read_id: String,
    /// The location of the node on the read.
    pub location: Wtdbg2ReadLocation,
}

/// A location on a read from wtdbg2.
pub struct Wtdbg2ReadLocation {
    /// The direction of the read location. True is forwards, false is reverse.
    pub direction: bool,
    /// The first bucket of the read location.
    pub bucket_offset: usize,
    /// The length in buckets of the read location.
    pub bucket_len: usize,
}

impl<'a> From<&'a str> for PlainWtdbg2NodeData {
    fn from(string: &'a str) -> Self {
        let mut split = string.split('\t');
        let id = split.next().unwrap()[1..].parse().unwrap();
        split.next();
        let mut read_associations = Vec::new();

        for read_association in split {
            let mut split = read_association.rsplitn(4, '_');
            let bucket_len = if let Ok(bucket_len) = split.next().unwrap().parse() {
                bucket_len
            } else {
                continue;
            };
            let bucket_offset = split.next().unwrap().parse().unwrap();
            let direction = split.next().unwrap();
            let direction = match direction {
                "F" => true,
                "R" => false,
                unknown => panic!("Unknown direction: '{}'", unknown),
            };
            let read_id = split.next().unwrap().to_owned();

            read_associations.push(Wtdbg2ReadAssociation {
                read_id,
                location: Wtdbg2ReadLocation {
                    direction,
                    bucket_offset,
                    bucket_len,
                },
            })
        }

        PlainWtdbg2NodeData {
            index: id,
            read_associations,
        }
    }
}

impl<'a> From<&'a str> for Wtdbg2ReadLocation {
    fn from(string: &'a str) -> Self {
        let mut split = string.split('_');
        let direction = split.next().unwrap();
        let direction = match direction {
            "F" => true,
            "R" => false,
            unknown => panic!("Unknown direction: '{}'", unknown),
        };
        let bucket_offset = split.next().unwrap().parse().unwrap();
        let bucket_len = split.next().unwrap().parse().unwrap();

        Self {
            direction,
            bucket_offset,
            bucket_len,
        }
    }
}

impl Wtdbg2ReadLocation {
    /// Returns true if the two locations on the same read form an edge in the fuzzy de Bruijn graph of wtdbg2.
    pub fn forms_edge(&self, head: &Self) -> bool {
        self.direction == head.direction
            && self.bucket_offset < head.bucket_offset
            && self.bucket_offset + self.bucket_len >= head.bucket_offset
    }
}

/// Read a genome graph in wtdbg2 format from a set of files.
pub fn read_graph_from_wtdbg2_from_files<
    P1: AsRef<Path>,
    P2: AsRef<Path>,
    NodeData: From<PlainWtdbg2NodeData>,
    EdgeData: Default + Clone,
    Graph: DynamicGraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    nodes_file: P1,
    reads_file: P2,
) -> Result<Graph> {
    read_graph_from_wtdbg2(
        BufReader::new(File::open(nodes_file)?),
        BufReader::new(File::open(reads_file)?),
    )
}

/// Read a genome graph in wtdbg2 format from a set of `BufRead`s.
pub fn read_graph_from_wtdbg2<
    R1: BufRead,
    R2: BufRead,
    NodeData: From<PlainWtdbg2NodeData>,
    EdgeData: Default + Clone,
    Graph: DynamicGraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    nodes: R1,
    reads: R2,
) -> Result<Graph> {
    let mut graph = Graph::default();

    info!("Loading nodes");
    for line in nodes.lines() {
        let line = line?;
        graph.add_node(PlainWtdbg2NodeData::from(line.as_str()).into());
    }
    info!("Loaded {} nodes", graph.node_count());

    info!("Loading edges");
    for line in reads.lines() {
        let line = line?;
        let mut split = line.split('\t');
        let _read_id = split.next().unwrap().to_owned();
        split.next();
        split.next();

        let mut last_node = split.next().unwrap();
        for current_node in split {
            let mut last_split = last_node.split(':');
            let mut current_split = current_node.split(':');

            let last_node_index: usize =
                if let Ok(last_node_id) = last_split.next().unwrap()[1..].parse() {
                    last_node_id
                } else {
                    continue;
                };
            let last_read_location = Wtdbg2ReadLocation::from(last_split.next().unwrap());
            let current_node_index: usize =
                if let Ok(current_node_id) = current_split.next().unwrap()[1..].parse() {
                    current_node_id
                } else {
                    continue;
                };
            let current_read_location = Wtdbg2ReadLocation::from(current_split.next().unwrap());

            if last_read_location.forms_edge(&current_read_location) {
                if graph.contains_edge_between(last_node_index.into(), current_node_index.into()) {
                    todo!("Edge already exists, increase count and add read association")
                } else {
                    todo!("Edge does not exist")
                }
            }

            last_node = current_node;
        }
    }
    info!("Loaded {} edges", graph.edge_count());

    todo!()
}

/// Write a list of contigs in wtdbg's .ctg.lay format.
pub fn write_contigs_to_wtdbg2<R: Read>(raw_reads: bio::io::fasta::Reader<R>) -> Result<()> {
    let mut read_map = HashMap::new();
    info!("Loading reads");
    for record in raw_reads.records() {
        let record = record?;
        read_map.insert(record.id().to_owned(), record.seq().to_owned());
    }
    info!("Loaded {} reads", read_map.len());
    todo!()
}
