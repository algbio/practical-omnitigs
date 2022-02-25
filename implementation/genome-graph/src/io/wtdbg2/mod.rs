use crate::error::Result;
use bigraph::interface::dynamic_bigraph::DynamicBigraph;
use bigraph::interface::BidirectedData;
use bigraph::traitgraph::interface::{Edge, ImmutableGraphContainer, StaticGraph};
use bigraph::traitgraph::traitsequence::interface::Sequence;
use bigraph::traitgraph::walks::{EdgeWalk, VecNodeWalk};
use compact_genome::implementation::DefaultGenome;
use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use compact_genome::interface::sequence::GenomeSequence;
use compact_genome::interface::sequence::{EditableGenomeSequence, OwnedGenomeSequence};
use regex::Regex;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt::{Debug, Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::time::{Duration, Instant};

/// Reading and writing the dot format of wtdbg2.
pub mod dot;

/// Node data as given in a .1.nodes file from wtdbg2.
#[derive(Clone, Debug)]
pub struct PlainWtdbg2NodeData {
    /// The index of the node in wtdbg2.
    pub index: usize,
    /// True if this is the variant of the node as given in the .1.nodes file, false if it is its reverse.
    pub forward: bool,
    /// True if the node is closed in wtdbg2.
    pub closed: bool,
    /// The read associations of the node.
    pub read_associations: Vec<Wtdbg2NodeReadAssociation>,
}

/// Read associations of nodes of wtdbg2.
#[derive(Clone, Debug)]
pub struct Wtdbg2NodeReadAssociation {
    /// The identifier of the read.
    pub read_id: String,
    /// The location of the node on the read.
    pub location: Wtdbg2ReadLocation,
    /// True if this read location contains a star.
    pub star: bool,
    /// True if this read location contains an exclamation mark.
    pub exclamation_mark: bool,
}

/// A location on a read from wtdbg2.
#[derive(Clone, Eq, PartialEq, Debug)]
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
        let id = &split.next().unwrap()[1..];
        let closed = id.ends_with('*');
        let id = if closed {
            id[..id.len() - 1].parse()
        } else {
            id.parse()
        }
        .unwrap();
        split.next();
        let mut read_associations = Vec::new();

        for read_association in split {
            let mut split = read_association.rsplitn(4, '_');
            let bucket_len = split.next().unwrap();
            let exclamation_mark = bucket_len.ends_with('!');
            let star = bucket_len.ends_with('*');
            let bucket_len = if exclamation_mark || star {
                bucket_len[..bucket_len.len() - 1].parse()
            } else {
                bucket_len.parse()
            }
            .unwrap();
            let bucket_offset = split.next().unwrap().parse().unwrap();
            let direction = split.next().unwrap();
            let direction = match direction {
                "F" => true,
                "R" => false,
                unknown => panic!("Unknown direction: '{}'", unknown),
            };
            let read_id = split.next().unwrap().to_owned();

            read_associations.push(Wtdbg2NodeReadAssociation {
                read_id,
                location: Wtdbg2ReadLocation {
                    direction,
                    bucket_offset,
                    bucket_len,
                },
                star,
                exclamation_mark,
            })
        }

        PlainWtdbg2NodeData {
            index: id,
            forward: true,
            closed,
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

impl PlainWtdbg2NodeData {
    /// Clone this node data for the reverse node.
    pub fn clone_reverse(&self) -> Self {
        Self {
            index: self.index,
            forward: !self.forward,
            closed: self.closed,
            read_associations: self
                .read_associations
                .iter()
                .map(|ra| {
                    let mut ra = ra.clone();
                    ra.location.direction = !ra.location.direction;
                    ra
                })
                .collect(),
        }
    }
}

impl Wtdbg2ReadLocation {
    /// Returns true if the two locations on the same read form an edge in the fuzzy de Bruijn graph of wtdbg2.
    pub fn forms_edge(&self, _head: &Self) -> bool {
        //self.direction == head.direction
        //    && self.bucket_offset < head.bucket_offset
        //    && self.bucket_offset + self.bucket_len >= head.bucket_offset
        // probably everything that is consecutive on a read forms an edge.
        true
    }
}

/// Functionalities of the wtdbg2 node data.
pub trait Wtdbg2NodeData {
    /// Returns the index of the node.
    fn index(&self) -> usize;

    /// Returns true if the node is forwards, and false if it is backwards.
    /// Forwards is always the node as given in the .1.nodes file.
    fn forward(&self) -> bool;
}

impl Wtdbg2NodeData for PlainWtdbg2NodeData {
    fn index(&self) -> usize {
        self.index
    }

    fn forward(&self) -> bool {
        self.forward
    }
}

/// Edge data derived from a .1.reads file.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct PlainWtdbg2EdgeData {
    /// The read associations of the edge. That are the locations on the reads that give evidence for this edge.
    pub read_associations: Vec<Wtdbg2EdgeReadAssociation>,
}

/// Read associations of edges of wtdbg2.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Wtdbg2EdgeReadAssociation {
    /// The identifier of the read.
    pub read_id: String,
    /// The location of the node on the read.
    pub location: Wtdbg2ReadLocation,
    /// True if the read location on the from node contains a star.
    pub from_star: bool,
    /// True if the read location on the from node contains an exclamation mark.
    pub from_exclamation_mark: bool,
    /// True if the read location on the to node contains a star.
    pub to_star: bool,
    /// True if the read location on the to node contains an exclamation mark.
    pub to_exclamation_mark: bool,
}

/// Functionalities of the wtdbg2 edge data.
/// Since when reading from wtdbg2 into a graph the edge data is converted into the target type before the edge is complete, the target type must implement this trait to allow adding to an edge after conversion into the target type.
pub trait Wtdbg2EdgeData {
    /// Add a new read association to existing edge data.
    fn add_edge_read_association(&mut self, read_association: Wtdbg2EdgeReadAssociation);

    /// Returns the read associations of the edge.
    fn edge_read_associations(&self) -> &[Wtdbg2EdgeReadAssociation];

    /// Returns the multiplicity of this edge, that is how many reads give evidence for this edge.
    fn multiplicity(&self) -> usize;

    /// Returns the length of this edge in buckets.
    /// That is typically the median of the lengths of the read fragments supporting the edge.
    fn length(&self) -> usize;
}

impl Wtdbg2EdgeData for PlainWtdbg2EdgeData {
    fn add_edge_read_association(&mut self, read_association: Wtdbg2EdgeReadAssociation) {
        self.read_associations.push(read_association);
    }

    fn edge_read_associations(&self) -> &[Wtdbg2EdgeReadAssociation] {
        &self.read_associations
    }

    fn multiplicity(&self) -> usize {
        self.read_associations.len()
    }

    fn length(&self) -> usize {
        let mut lengths = Vec::new();
        for read_association in &self.read_associations {
            lengths.push(read_association.location.bucket_len);
        }
        lengths.sort_unstable();
        lengths[(lengths.len().max(1) - 1) / 2]
    }
}

impl BidirectedData for PlainWtdbg2EdgeData {
    fn mirror(&self) -> Self {
        Self {
            read_associations: self
                .read_associations
                .iter()
                .map(|ra| Wtdbg2EdgeReadAssociation {
                    read_id: ra.read_id.clone(),
                    location: ra.location.mirror(),
                    from_star: ra.to_star,
                    from_exclamation_mark: ra.to_exclamation_mark,
                    to_star: ra.from_star,
                    to_exclamation_mark: ra.from_exclamation_mark,
                })
                .collect(),
        }
    }
}

impl BidirectedData for Wtdbg2ReadLocation {
    fn mirror(&self) -> Self {
        Self {
            direction: !self.direction,
            bucket_offset: self.bucket_offset,
            bucket_len: self.bucket_len,
        }
    }
}

/// Read a genome graph in wtdbg2 format from a set of files.
pub fn read_graph_from_wtdbg2_from_files<
    P1: AsRef<Path>,
    P2: AsRef<Path>,
    P3: AsRef<Path>,
    NodeData: From<PlainWtdbg2NodeData> + Wtdbg2NodeData,
    EdgeData: From<PlainWtdbg2EdgeData> + Wtdbg2EdgeData,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    nodes_file: P1,
    reads_file: P2,
    dot_file: P3,
) -> Result<Graph> {
    read_graph_from_wtdbg2(
        BufReader::new(File::open(nodes_file)?),
        BufReader::new(File::open(reads_file)?),
        BufReader::new(File::open(dot_file)?),
    )
}

/// Read a genome graph in wtdbg2 format from a set of `BufRead`s.
pub fn read_graph_from_wtdbg2<
    R1: BufRead,
    R2: BufRead,
    R3: BufRead,
    NodeData: From<PlainWtdbg2NodeData> + Wtdbg2NodeData,
    EdgeData: From<PlainWtdbg2EdgeData> + Wtdbg2EdgeData,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    nodes: R1,
    reads: R2,
    dot: R3,
) -> Result<Graph> {
    let mut graph = Graph::default();
    let mut node_map = HashMap::new();

    info!("Loading nodes");
    for line in nodes.lines() {
        let line = line?;
        let forward_node_data = PlainWtdbg2NodeData::from(line.as_str());
        let reverse_node = graph.add_node(forward_node_data.clone_reverse().into());
        let forward_node = graph.add_node(forward_node_data.clone().into());
        node_map.insert(forward_node_data.index, (reverse_node, forward_node));
        graph.set_mirror_nodes(forward_node, reverse_node);
    }
    info!("Loaded {} nodes", graph.node_count());

    info!("Loading edges from dot");
    for line in dot.lines() {
        let line = line?;
        if !line.contains("->") {
            continue;
        }

        let mut split = line.split(' ');
        let n1 = node_map
            .get(&split.next().unwrap()[1..].parse().unwrap())
            .unwrap();
        split.next();
        let n2 = node_map
            .get(&split.next().unwrap()[1..].parse().unwrap())
            .unwrap();
        let label = &split.next().unwrap()[8..10];
        let from_forward = match label.as_bytes()[0] {
            b'+' => true,
            b'-' => false,
            unknown => bail!("Unknown node direction: {}", unknown),
        };
        let to_forward = match label.as_bytes()[1] {
            b'+' => true,
            b'-' => false,
            unknown => bail!("Unknown node direction: {}", unknown),
        };

        let n1 = if from_forward { n1.1 } else { n1.0 };
        let n2 = if to_forward { n2.1 } else { n2.0 };
        graph.add_edge(
            n1,
            n2,
            PlainWtdbg2EdgeData {
                read_associations: Vec::new(),
            }
            .into(),
        );
    }

    info!("Loading edge read associations");
    for line in reads.lines() {
        let line = line?;
        let mut split = line.split('\t');
        let read_id = split.next().unwrap().to_owned();
        split.next();
        let node_association_amount: usize = split.next().unwrap().parse().unwrap();
        if node_association_amount < 2 {
            continue;
        }
        let nodes: Vec<_> = split.collect();

        for (index, n1) in nodes.iter().enumerate() {
            let mut n1_split = n1.split(':');
            let n1_index = &n1_split.next().unwrap()[1..];
            let n1_star = n1_index.ends_with('*');
            let n1_exclamation_mark = n1_index.ends_with('!');
            let n1_index = node_map
                .get(
                    &if n1_star || n1_exclamation_mark {
                        n1_index[..n1_index.len() - 1].parse()
                    } else {
                        n1_index.parse()
                    }
                    .unwrap(),
                )
                .unwrap();
            let n1_read_location = Wtdbg2ReadLocation::from(n1_split.next().unwrap());

            let n1_node_index = if n1_read_location.direction {
                n1_index.1
            } else {
                n1_index.0
            };
            let reverse_n2_node_index = if n1_read_location.direction {
                n1_index.0
            } else {
                n1_index.1
            };

            for n2 in nodes.iter().skip(index + 1) {
                let mut n2_split = n2.split(':');
                let n2_index = &n2_split.next().unwrap()[1..];
                let n2_star = n2_index.ends_with('*');
                let n2_exclamation_mark = n2_index.ends_with('!');
                let n2_index = node_map
                    .get(
                        &if n2_star || n2_exclamation_mark {
                            n2_index[..n2_index.len() - 1].parse()
                        } else {
                            n2_index.parse()
                        }
                        .unwrap(),
                    )
                    .unwrap();
                let n2_read_location = Wtdbg2ReadLocation::from(n2_split.next().unwrap());

                let n2_node_index = if n2_read_location.direction {
                    n2_index.1
                } else {
                    n2_index.0
                };
                let reverse_n1_node_index = if n2_read_location.direction {
                    n2_index.0
                } else {
                    n2_index.1
                };

                let existing_edge = graph.edges_between(n1_node_index, n2_node_index).next();
                if let Some(edge) = existing_edge {
                    let forward_read_association = Wtdbg2EdgeReadAssociation {
                        read_id: read_id.clone(),
                        location: Wtdbg2ReadLocation {
                            direction: true,
                            bucket_offset: n1_read_location.bucket_offset,
                            bucket_len: n2_read_location.bucket_offset
                                + n2_read_location.bucket_len
                                - n1_read_location.bucket_offset,
                        },
                        from_star: n1_star,
                        to_star: n2_star,
                        from_exclamation_mark: n1_exclamation_mark,
                        to_exclamation_mark: n2_exclamation_mark,
                    };

                    let edge_data = graph.edge_data_mut(edge);
                    edge_data.add_edge_read_association(forward_read_association);
                }

                let existing_edge = graph
                    .edges_between(reverse_n1_node_index, reverse_n2_node_index)
                    .next();
                if let Some(edge) = existing_edge {
                    let reverse_read_association = Wtdbg2EdgeReadAssociation {
                        read_id: read_id.clone(),
                        location: Wtdbg2ReadLocation {
                            direction: false,
                            bucket_offset: n1_read_location.bucket_offset,
                            bucket_len: n2_read_location.bucket_offset
                                + n2_read_location.bucket_len
                                - n1_read_location.bucket_offset,
                        },
                        from_star: n2_star,
                        to_star: n1_star,
                        from_exclamation_mark: n2_exclamation_mark,
                        to_exclamation_mark: n1_exclamation_mark,
                    };

                    let edge_data = graph.edge_data_mut(edge);
                    edge_data.add_edge_read_association(reverse_read_association);
                }
            }
        }
    }

    info!("Loaded {} edges", graph.edge_count());

    for edge_index in graph.edge_indices() {
        let edge_data = graph.edge_data(edge_index);
        if edge_data.edge_read_associations().is_empty() {
            let from_node = graph.edge_endpoints(edge_index).from_node;
            let to_node = graph.edge_endpoints(edge_index).to_node;
            warn!(
                "Edge has no read associations: N{}{} -> N{}{}",
                graph.node_data(from_node).index(),
                if graph.node_data(from_node).forward() {
                    '+'
                } else {
                    '-'
                },
                graph.node_data(to_node).index(),
                if graph.node_data(to_node).forward() {
                    '+'
                } else {
                    '-'
                },
            );
        } else if edge_data.edge_read_associations().len() < 3 {
            let from_node = graph.edge_endpoints(edge_index).from_node;
            let to_node = graph.edge_endpoints(edge_index).to_node;
            debug!(
                "Edge has only {} read associations: N{}{} -> N{}{}",
                edge_data.edge_read_associations().len(),
                graph.node_data(from_node).index(),
                if graph.node_data(from_node).forward() {
                    '+'
                } else {
                    '-'
                },
                graph.node_data(to_node).index(),
                if graph.node_data(to_node).forward() {
                    '+'
                } else {
                    '-'
                },
            );
        }
    }

    Ok(graph)
}

/// A .ctg.lay file.
#[derive(Default, Debug, Clone)]
pub struct RawWtdbg2Contigs {
    contigs: Vec<RawWtdbg2Contig>,
}

/// A contig from a .ctg.lay file.
#[derive(Debug, Clone)]
pub struct RawWtdbg2Contig {
    index: usize,
    nodes: usize,
    len: usize,
    edges: Vec<RawWtdbg2ContigEdge>,
}

/// An edge from a .ctg.lay file.
#[derive(Debug, Clone)]
pub struct RawWtdbg2ContigEdge {
    offset: usize,
    from_node: usize,
    from_direction: bool,
    to_node: usize,
    to_direction: bool,
    supports: Vec<RawWtdbg2ContigEdgeSupport>,
}

/// An edge support from a .ctg.lay file.
#[derive(Clone, Debug)]
pub struct RawWtdbg2ContigEdgeSupport {
    read: String,
    direction: bool,
    offset: usize,
    len: usize,
    sequence: String,
}

impl RawWtdbg2Contigs {
    /// Sets the indices of the contigs according to their order.
    pub fn update_indices(&mut self) {
        for (index, contig) in self.contigs.iter_mut().enumerate() {
            contig.index = index + 1
        }
    }

    /// Sorts the contigs lexicographically by their node ids, and sorts the edge supports uniquely as well.
    pub fn sort_contigs_topologically(&mut self) {
        self.contigs.sort_unstable_by(|a, b| {
            for (ea, eb) in a.edges.iter().zip(b.edges.iter()) {
                if ea.from_node != eb.from_node {
                    return ea.from_node.cmp(&eb.from_node);
                }
            }

            if let Some((ea, eb)) = a.edges.iter().zip(b.edges.iter()).last() {
                if ea.to_node != eb.to_node {
                    return ea.to_node.cmp(&eb.to_node);
                }
            }

            if a.edges.len() != b.edges.len() {
                return a.edges.len().cmp(&b.edges.len());
            }

            Ordering::Equal
        });

        for contig in self.contigs.iter_mut() {
            for edge in contig.edges.iter_mut() {
                edge.supports.sort_unstable_by(|a, b| {
                    (&a.read, a.direction, a.len, &a.sequence).cmp(&(
                        &b.read,
                        b.direction,
                        b.len,
                        &b.sequence,
                    ))
                });
            }
        }
    }

    /// Print differences between the two contig lists to stdout.
    pub fn compare_contigs(&self, other: &Self) {
        if self.contigs.len() != other.contigs.len() {
            warn!(
                "different number of contigs: {} != {}",
                self.contigs.len(),
                other.contigs.len()
            );
        }

        for (a, b) in self.contigs.iter().zip(other.contigs.iter()) {
            if a.index != b.index {
                warn!("contigs differ in index: {} != {}", a.index, b.index);
            }
            if a.len != b.len {
                warn!("ctg{} differs in length: {} != {}", a.index, a.len, b.len);
            }
            if a.nodes != b.nodes {
                error!(
                    "ctg{} differs in nodes: {} != {}",
                    a.index, a.nodes, b.nodes
                );
            }
            if a.edges.len() != b.edges.len() {
                error!(
                    "ctg{} differs in edge length: {} != {}",
                    a.index,
                    a.edges.len(),
                    b.edges.len()
                );
            }

            let mut node_diff = false;
            for (i, (ea, eb)) in a.edges.iter().zip(b.edges.iter()).enumerate() {
                if node_diff {
                    continue;
                }
                if ea.offset != eb.offset {
                    warn!(
                        "ctg{} edge {} differs in offset: {} != {}",
                        a.index, i, ea.offset, eb.offset
                    );
                }
                if ea.from_node != eb.from_node {
                    error!(
                        "ctg{} edge {} differs in from node: {} != {}",
                        a.index, i, ea.from_node, eb.from_node
                    );
                    node_diff = true;
                }
                if ea.from_direction != eb.from_direction {
                    error!(
                        "ctg{} edge {} differs in from direction: {} != {}",
                        a.index, i, ea.from_direction, eb.from_direction
                    );
                    node_diff = true;
                }
                if ea.to_node != eb.to_node {
                    error!(
                        "ctg{} edge {} differs in to node: {} != {}",
                        a.index, i, ea.to_node, eb.to_node
                    );
                    node_diff = true;
                }
                if ea.to_direction != eb.to_direction {
                    error!(
                        "ctg{} edge {} differs in from node: {} != {}",
                        a.index, i, ea.to_direction, eb.to_direction
                    );
                    node_diff = true;
                }
                if node_diff {
                    continue;
                }
                if ea.supports.len() != eb.supports.len() {
                    error!(
                        "ctg{} edge {} differs in support amount: {} != {}",
                        a.index,
                        i,
                        ea.supports.len(),
                        eb.supports.len()
                    );
                }
                for (si, (sa, sb)) in ea.supports.iter().zip(eb.supports.iter()).enumerate() {
                    let mut read_dir_diff = false;
                    if sa.read != sb.read {
                        error!(
                            "ctg{} edge {} support {} differs in read: {} != {}",
                            a.index, i, si, sa.read, sb.read
                        );
                        read_dir_diff = true;
                    }
                    if sa.direction != sb.direction {
                        error!(
                            "ctg{} edge {} support {} differs in direction: {} != {}",
                            a.index, i, si, sa.direction, sb.direction
                        );
                        read_dir_diff = true;
                    }
                    if sa.offset != sb.offset {
                        error!(
                            "ctg{} edge {} support {} differs in offset: {} != {}",
                            a.index, i, si, sa.offset, sb.offset
                        );
                    }
                    if sa.len != sb.len {
                        error!(
                            "ctg{} edge {} support {} differs in length: {} != {}",
                            a.index, i, si, sa.len, sb.len
                        );
                    }
                    if sa.sequence != sb.sequence && !read_dir_diff {
                        error!(
                            "ctg{} edge {} support {} differs in sequence: {} != {}",
                            a.index, i, si, sa.sequence, sb.sequence
                        );
                    }
                }
            }
        }
    }
}

impl Display for RawWtdbg2Contigs {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for contig in &self.contigs {
            Display::fmt(contig, f)?;
        }
        Ok(())
    }
}

impl Display for RawWtdbg2Contig {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "<ctg{} nodes={} len={}",
            self.index, self.nodes, self.len
        )?;
        for edge in &self.edges {
            Display::fmt(edge, f)?;
        }
        Ok(())
    }
}

impl Display for RawWtdbg2ContigEdge {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "E\t{}\tN{}\t{}\tN{}\t{}",
            self.offset,
            self.from_node,
            if self.from_direction { '+' } else { '-' },
            self.to_node,
            if self.to_direction { '+' } else { '-' }
        )?;
        for support in &self.supports {
            Display::fmt(support, f)?;
        }
        Ok(())
    }
}

impl Display for RawWtdbg2ContigEdgeSupport {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "S\t{}\t{}\t{}\t{}\t{}",
            self.read,
            if self.direction { '+' } else { '-' },
            self.offset,
            self.len,
            self.sequence
        )
    }
}

/// Convert a list of walks into a RawWtdbg2Contigs struct that represents a .ctg.lay file.
/// This opens the given path as raw reads file in fasta format.
pub fn convert_walks_to_wtdbg2_contigs_with_file<
    'ws,
    P: AsRef<Path> + Debug,
    NodeData: Wtdbg2NodeData,
    EdgeData: Wtdbg2EdgeData,
    Graph: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
>(
    graph: &Graph,
    walks: WalkSource,
    raw_reads_file: P,
) -> Result<RawWtdbg2Contigs> {
    convert_walks_to_wtdbg2_contigs(
        graph,
        walks,
        bio::io::fasta::Reader::from_file(raw_reads_file)?,
    )
}

/// Convert a list of walks into a RawWtdbg2Contigs struct that represents a .ctg.lay file.
/// This interprets the given reader as raw reads source in fasta format.
pub fn convert_walks_to_wtdbg2_contigs<
    'ws,
    R: BufRead,
    NodeData: Wtdbg2NodeData,
    EdgeData: Wtdbg2EdgeData,
    Graph: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
>(
    graph: &Graph,
    walks: WalkSource,
    raw_reads: bio::io::fasta::Reader<R>,
) -> Result<RawWtdbg2Contigs> {
    let mut read_map = HashMap::<_, DefaultGenome<DnaAlphabet>>::new();
    info!("Loading reads");
    let mut last_print_time = Instant::now();

    let trim_regex = Regex::new(r"^(.+?)/[0-9]+_[0-9]+$").unwrap();
    for record in raw_reads.records() {
        if Instant::now() - last_print_time >= Duration::from_secs(5) {
            info!("Loaded {} reads", read_map.len());
            last_print_time = Instant::now();
        }

        let record = record?;
        // For some reason, wtdbg2 trims the last parts of reads.
        // Maybe these indicate positions on an actual read.
        // If wtdbg2 combines the reads then this is wrong.
        let id = if trim_regex.is_match(record.id()) {
            &record.id()[..record.id().rfind('/').unwrap()]
        } else {
            record.id()
        }
        .to_owned();
        if let Some(genome) = read_map.get_mut(&id) {
            genome.extend_from_slice_u8(record.seq()).unwrap();
        } else {
            read_map.insert(
                id.to_owned(),
                DefaultGenome::from_slice_u8(record.seq()).unwrap(),
            );
        }
    }
    info!("Finished loading {} reads", read_map.len());

    let walk_iter = walks.into_iter();
    if let (min, Some(max)) = walk_iter.size_hint() {
        if min == max {
            info!("Converting {} walks", min);
        } else {
            info!("Converting between {} and {} walks", min, max);
        }
    } else {
        info!("Converting at least {} walks", walk_iter.size_hint().0);
    }
    last_print_time = Instant::now();
    let mut dropped_walks = 0usize;
    let mut printed_walks = 0usize;
    let mut contigs = RawWtdbg2Contigs::default();

    for walk in walk_iter {
        let walk_index = printed_walks;
        if Instant::now() - last_print_time >= Duration::from_secs(5) {
            info!("Converted {} walks", walk_index);
            last_print_time = Instant::now();
        }

        // Compute edge offsets on the walk.
        // The length of an edge is the median length of its supporting read fragments.
        let mut offsets = vec![0];

        for &edge in walk.iter() {
            offsets.push(offsets.last().unwrap() + graph.edge_data(edge).length() - 4);
        }

        if offsets.last().unwrap() * 256 < 5000 || walk.len() < 2 {
            dropped_walks += 1;
            continue;
        } else {
            printed_walks += 1;
        }

        let mut edges = Vec::new();

        for (&edge, offset) in walk.iter().zip(offsets.iter()) {
            let Edge { from_node, to_node } = graph.edge_endpoints(edge);
            let from_node_data = graph.node_data(from_node);
            let to_node_data = graph.node_data(to_node);

            let mut read_associations = Vec::new();
            for read_association in graph.edge_data(edge).edge_read_associations() {
                let offset = read_association.location.bucket_offset * 256;
                let len = read_association.location.bucket_len * 256;
                read_associations.push(RawWtdbg2ContigEdgeSupport {
                    read: read_association.read_id.to_owned(),
                    direction: read_association.location.direction,
                    offset,
                    len,
                    sequence: read_map.get(&read_association.read_id).unwrap()
                        [offset..offset + len]
                        .as_string(),
                });
            }

            edges.push(RawWtdbg2ContigEdge {
                offset: offset * 256,
                from_node: from_node_data.index(),
                from_direction: from_node_data.forward(),
                to_node: to_node_data.index(),
                to_direction: to_node_data.forward(),
                supports: read_associations,
            });
        }

        contigs.contigs.push(RawWtdbg2Contig {
            index: walk_index + 1,
            nodes: walk.len() + 1,
            len: offsets.last().unwrap() * 256,
            edges,
        });
    }

    info!("Finished converting {} walks", printed_walks);
    info!("{} too short walks were dropped", dropped_walks);

    Ok(contigs)
}

/// Read a .ctg.lay file into a RawWtdbg2Contigs struct.
pub fn read_wtdbg2_contigs_from_file<P: AsRef<Path>>(input_file: P) -> Result<RawWtdbg2Contigs> {
    read_wtdbg2_contigs(File::open(input_file)?)
}

/// Read a .ctg.lay source into a RawWtdbg2Contigs struct.
pub fn read_wtdbg2_contigs<R: Read>(input: R) -> Result<RawWtdbg2Contigs> {
    let mut result = RawWtdbg2Contigs::default();

    for line in BufReader::new(input).lines() {
        let line = line.unwrap();
        match line.as_bytes()[0] {
            b'>' => {
                let mut split = line.split(' ');
                let index = split.next().unwrap()[4..].parse().unwrap();
                let nodes = split.next().unwrap()[6..].parse().unwrap();
                let len = split.next().unwrap()[4..].parse().unwrap();
                result.contigs.push(RawWtdbg2Contig {
                    index,
                    nodes,
                    len,
                    edges: Vec::new(),
                });
            }
            b'E' => {
                let mut split = line.split('\t');
                split.next().unwrap();
                let offset = split.next().unwrap().parse().unwrap();
                let from_node = split.next().unwrap()[1..].parse().unwrap();
                let from_direction = match split.next().unwrap() {
                    "+" => true,
                    "-" => false,
                    unknown => bail!("Unknown direction: '{}'", unknown),
                };
                let to_node = split.next().unwrap()[1..].parse().unwrap();
                let to_direction = match split.next().unwrap() {
                    "+" => true,
                    "-" => false,
                    unknown => bail!("Unknown direction: '{}'", unknown),
                };
                result
                    .contigs
                    .last_mut()
                    .unwrap()
                    .edges
                    .push(RawWtdbg2ContigEdge {
                        offset,
                        from_node,
                        from_direction,
                        to_node,
                        to_direction,
                        supports: Vec::new(),
                    });
            }
            b'S' => {
                let mut split = line.split('\t');
                split.next().unwrap();
                let read = split.next().unwrap().to_owned();
                let direction = match split.next().unwrap() {
                    "+" => true,
                    "-" => false,
                    unknown => bail!("Unknown direction: '{}'", unknown),
                };
                let offset = split.next().unwrap().parse().unwrap();
                let len = split.next().unwrap().parse().unwrap();
                let sequence = split.next().unwrap().to_owned();

                result
                    .contigs
                    .last_mut()
                    .unwrap()
                    .edges
                    .last_mut()
                    .unwrap()
                    .supports
                    .push(RawWtdbg2ContigEdgeSupport {
                        read,
                        direction,
                        offset,
                        len,
                        sequence,
                    });
            }
            error => bail!("Unknown line start: '{}'", error as char),
        }
    }

    Ok(result)
}

/// Write a list of contigs in wtdbg's .ctg.lay format to a file.
pub fn write_contigs_to_wtdbg2_to_file<
    'ws,
    P1: AsRef<Path> + Debug,
    P2: AsRef<Path>,
    NodeData: Wtdbg2NodeData,
    EdgeData: Wtdbg2EdgeData,
    Graph: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
>(
    graph: &Graph,
    walks: WalkSource,
    raw_reads_file: P1,
    output_file: P2,
) -> Result<()> {
    write_contigs_to_wtdbg2(
        graph,
        walks,
        bio::io::fasta::Reader::from_file(raw_reads_file)?,
        &mut BufWriter::new(File::create(output_file)?),
    )
}

/// Write a list of contigs in wtdbg's .ctg.lay format.
pub fn write_contigs_to_wtdbg2<
    'ws,
    R: BufRead,
    W: Write,
    NodeData: Wtdbg2NodeData,
    EdgeData: Wtdbg2EdgeData,
    Graph: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
>(
    graph: &Graph,
    walks: WalkSource,
    raw_reads: bio::io::fasta::Reader<R>,
    output: &mut W,
) -> Result<()> {
    let mut read_map = HashMap::<_, DefaultGenome<DnaAlphabet>>::new();
    info!("Loading reads");
    let mut last_print_time = Instant::now();

    let trim_regex = Regex::new(r"^(.+?)/[0-9]+_[0-9]+$").unwrap();
    for record in raw_reads.records() {
        if Instant::now() - last_print_time >= Duration::from_secs(5) {
            info!("Loaded {} reads", read_map.len());
            last_print_time = Instant::now();
        }

        let record = record?;
        // For some reason, wtdbg2 trims the last parts of reads.
        // Maybe these indicate positions on an actual read.
        // If wtdbg2 combines the reads then this is wrong.
        let id = if trim_regex.is_match(record.id()) {
            &record.id()[..record.id().rfind('/').unwrap()]
        } else {
            record.id()
        }
        .to_owned();
        if let Some(genome) = read_map.get_mut(&id) {
            genome.extend_from_slice_u8(record.seq()).unwrap();
        } else {
            read_map.insert(
                id.to_owned(),
                DefaultGenome::from_slice_u8(record.seq()).unwrap(),
            );
        }
    }
    info!("Finished loading {} reads", read_map.len());

    let mut walks: Vec<_> = walks
        .into_iter()
        .map(|w| {
            let mut len = 0;
            for &edge in w.iter() {
                len += graph.edge_data(edge).length();
            }
            len -= 4 * (w.len() - 1);
            (len, w)
        })
        .collect();
    walks.sort_by(|(l1, _), (l2, _)| l2.cmp(l1));

    let walk_iter = walks.iter().map(|(_, w)| w);
    if let (min, Some(max)) = walk_iter.size_hint() {
        if min == max {
            info!("Writing {} walks", min);
        } else {
            info!("Writing between {} and {} walks", min, max);
        }
    } else {
        info!("Writing at least {} walks", walk_iter.size_hint().0);
    }
    last_print_time = Instant::now();
    let mut dropped_walks = 0usize;
    let mut printed_walks = 0usize;

    for walk in walk_iter {
        let walk_index = printed_walks;
        if Instant::now() - last_print_time >= Duration::from_secs(5) {
            info!("Wrote {} walks", walk_index);
            last_print_time = Instant::now();
        }

        // Compute edge offsets on the walk.
        // The length of an edge is the median length of its supporting read fragments.
        let mut offsets = vec![0];

        for &edge in walk.iter() {
            offsets.push(offsets.last().unwrap() + graph.edge_data(edge).length() - 4);
        }

        if offsets.last().unwrap() * 256 < 5000 || walk.len() < 2 {
            dropped_walks += 1;
            continue;
        } else {
            printed_walks += 1;
        }

        writeln!(
            output,
            ">ctg{} nodes={} len={}",
            walk_index + 1,
            walk.len() + 1,
            offsets.last().unwrap() * 256
        )?;

        for (&edge, offset) in walk.iter().zip(offsets.iter()) {
            let Edge { from_node, to_node } = graph.edge_endpoints(edge);
            let from_node_data = graph.node_data(from_node);
            let to_node_data = graph.node_data(to_node);
            writeln!(
                output,
                "E\t{}\tN{}\t{}\tN{}\t{}",
                offset * 256,
                from_node_data.index(),
                if from_node_data.forward() { '+' } else { '-' },
                to_node_data.index(),
                if to_node_data.forward() { '+' } else { '-' }
            )?;

            for read_association in graph.edge_data(edge).edge_read_associations() {
                let offset = read_association.location.bucket_offset * 256;
                let len = read_association.location.bucket_len * 256;
                writeln!(
                    output,
                    "S\t{}\t{}\t{}\t{}\t{:?}",
                    read_association.read_id,
                    if read_association.location.direction {
                        '+'
                    } else {
                        '-'
                    },
                    offset,
                    len,
                    &read_map.get(&read_association.read_id).unwrap()[offset..offset + len]
                )?;
            }
        }
    }

    info!("Finished writing {} walks", printed_walks);
    info!("{} too short walks were dropped", dropped_walks);

    output.flush()?;
    Ok(())
}

/// Write a list of contigs as lists of wtdbg2's node ids to a file.
pub fn write_contigs_as_wtdbg2_node_ids_to_file<
    'ws,
    P: AsRef<Path>,
    NodeData: Wtdbg2NodeData,
    EdgeData: Wtdbg2EdgeData,
    Graph: StaticGraph<NodeData = NodeData, EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
>(
    graph: &Graph,
    walks: WalkSource,
    output_file: P,
) -> Result<()> {
    write_contigs_as_wtdbg2_node_ids(
        graph,
        walks,
        &mut BufWriter::new(File::create(output_file)?),
    )
}

/// Write a list of contigs as lists of wtdbg2's node ids.
pub fn write_contigs_as_wtdbg2_node_ids<
    'ws,
    W: Write,
    NodeData: Wtdbg2NodeData,
    EdgeData: Wtdbg2EdgeData,
    Graph: StaticGraph<NodeData = NodeData, EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph, Subwalk>,
    Subwalk: for<'w> EdgeWalk<'w, Graph, Subwalk> + ?Sized,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
>(
    graph: &Graph,
    walks: WalkSource,
    output: &mut W,
) -> Result<()> {
    for walk in walks {
        let walk: VecNodeWalk<Graph> = walk.clone_as_node_walk(graph).unwrap();
        for &node in walk.iter() {
            write!(
                output,
                "{} {} ",
                graph.node_data(node).index(),
                if graph.node_data(node).forward() {
                    "+"
                } else {
                    "-"
                }
            )?;
        }
        writeln!(output)?;
    }

    Ok(())
}

/// Build a graph from the unitigs of wtdbg2.
pub fn build_wtdbg2_unitigs_graph<
    Graph: DynamicBigraph<NodeData = (usize, bool), EdgeData = (RawWtdbg2Contig, bool)> + Default,
>(
    contigs: &RawWtdbg2Contigs,
) -> Graph {
    let mut graph = Graph::default();
    let mut node_map = HashMap::new();

    for contig in &contigs.contigs {
        let n1 = contig.edges.first().unwrap().from_node;
        let n1_dir = contig.edges.first().unwrap().from_direction;
        let n2 = contig.edges.first().unwrap().to_node;
        let n2_dir = contig.edges.first().unwrap().to_direction;

        let n1_index = if let Some(&index) = node_map.get(&(n1, n1_dir)) {
            index
        } else {
            let index = graph.add_node((n1, n1_dir));
            node_map.insert((n1, n1_dir), index);
            index
        };
        let n1_r_index = if let Some(&index) = node_map.get(&(n1, !n1_dir)) {
            index
        } else {
            let index = graph.add_node((n1, !n1_dir));
            node_map.insert((n1, !n1_dir), index);
            index
        };
        let n2_index = if let Some(&index) = node_map.get(&(n2, n2_dir)) {
            index
        } else {
            let index = graph.add_node((n2, n2_dir));
            node_map.insert((n2, n2_dir), index);
            index
        };
        let n2_r_index = if let Some(&index) = node_map.get(&(n2, !n2_dir)) {
            index
        } else {
            let index = graph.add_node((n2, !n2_dir));
            node_map.insert((n2, !n2_dir), index);
            index
        };

        graph.add_edge(n1_index, n2_index, (contig.clone(), true));
        graph.add_edge(n2_r_index, n1_r_index, (contig.clone(), false));
    }

    graph
}
