use crate::error::Result;
use bigraph::interface::dynamic_bigraph::DynamicBigraph;
use bigraph::interface::BidirectedData;
use bigraph::traitgraph::index::GraphIndex;
use bigraph::traitgraph::interface::{Edge, ImmutableGraphContainer};
use bigraph::traitgraph::walks::EdgeWalk;
use compact_genome::implementation::vector_genome_impl::VectorGenome;
use compact_genome::interface::ExtendableGenome;
use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::iter::FromIterator;
use std::path::Path;
use std::time::{Duration, Instant};

/// Node data as given in a .1.nodes file from wtdbg2.
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
#[derive(Clone)]
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
    NodeData: From<PlainWtdbg2NodeData>,
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
    NodeData: From<PlainWtdbg2NodeData>,
    EdgeData: From<PlainWtdbg2EdgeData> + Wtdbg2EdgeData,
    Graph: DynamicBigraph<NodeData = NodeData, EdgeData = EdgeData> + Default,
>(
    nodes: R1,
    reads: R2,
    dot: R3,
) -> Result<Graph> {
    let mut graph = Graph::default();

    info!("Loading nodes");
    for line in nodes.lines() {
        let line = line?;
        let forward_node_data = PlainWtdbg2NodeData::from(line.as_str());
        let reverse_node = graph.add_node(forward_node_data.clone_reverse().into());
        let forward_node = graph.add_node(forward_node_data.into());
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
        let n1: usize = split.next().unwrap()[1..].parse().unwrap();
        split.next();
        let n2: usize = split.next().unwrap()[1..].parse().unwrap();
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

        let n1 = 2 * n1 + if from_forward { 1 } else { 0 };
        let n2 = 2 * n2 + if to_forward { 1 } else { 0 };
        graph.add_edge(
            n1.into(),
            n2.into(),
            PlainWtdbg2EdgeData {
                read_associations: Vec::new(),
            }
            .into(),
        );
    }

    info!("Loading edges");
    for line in reads.lines() {
        let line = line?;
        let mut split = line.split('\t');
        let read_id = split.next().unwrap().to_owned();
        split.next();
        let node_association_amount: usize = split.next().unwrap().parse().unwrap();
        if node_association_amount < 2 {
            continue;
        }

        let mut last_node = split.next().unwrap();
        for current_node in split {
            let mut last_split = last_node.split(':');
            let mut current_split = current_node.split(':');

            let last_node_index = &last_split.next().unwrap()[1..];
            let last_star = last_node_index.ends_with('*');
            let last_exclamation_mark = last_node_index.ends_with('!');
            let last_node_index: usize = if last_star || last_exclamation_mark {
                last_node_index[..last_node_index.len() - 1].parse()
            } else {
                last_node_index.parse()
            }
            .unwrap();
            let last_read_location = Wtdbg2ReadLocation::from(last_split.next().unwrap());
            let current_node_index = &current_split.next().unwrap()[1..];
            let current_star = current_node_index.ends_with('*');
            let current_exclamation_mark = current_node_index.ends_with('!');
            let current_node_index: usize = if current_star || current_exclamation_mark {
                current_node_index[..current_node_index.len() - 1].parse()
            } else {
                current_node_index.parse()
            }
            .unwrap();
            let current_read_location = Wtdbg2ReadLocation::from(current_split.next().unwrap());

            if current_star || current_exclamation_mark {
                continue;
            }

            let from_node_index =
                (2 * last_node_index + if last_read_location.direction { 1 } else { 0 }).into();
            let to_node_index = (2 * current_node_index
                + if current_read_location.direction {
                    1
                } else {
                    0
                })
            .into();
            let reverse_from_node_index = (2 * current_node_index
                + if current_read_location.direction {
                    0
                } else {
                    1
                })
            .into();
            let reverse_to_node_index =
                (2 * last_node_index + if last_read_location.direction { 0 } else { 1 }).into();

            if last_read_location.forms_edge(&current_read_location) {
                let forward_read_association = Wtdbg2EdgeReadAssociation {
                    read_id: read_id.clone(),
                    location: Wtdbg2ReadLocation {
                        direction: true,
                        bucket_offset: last_read_location.bucket_offset,
                        bucket_len: current_read_location.bucket_offset
                            + current_read_location.bucket_len
                            - last_read_location.bucket_offset,
                    },
                    from_star: last_star,
                    to_star: current_star,
                    from_exclamation_mark: last_exclamation_mark,
                    to_exclamation_mark: current_exclamation_mark,
                };

                let existing_edge = graph.edges_between(from_node_index, to_node_index).next();
                if let Some(edge) = existing_edge {
                    let edge_data = graph.edge_data_mut(edge);
                    edge_data.add_edge_read_association(forward_read_association);
                } /*else {
                      let edge_data = PlainWtdbg2EdgeData {
                          read_associations: vec![forward_read_association],
                      };
                      graph.add_edge(from_node_index, to_node_index, edge_data.into());
                      warn!("Added non-existing edge N{} -> N{}", from_node_index.as_usize(), to_node_index.as_usize());
                  }*/

                let reverse_read_association = Wtdbg2EdgeReadAssociation {
                    read_id: read_id.clone(),
                    location: Wtdbg2ReadLocation {
                        direction: false,
                        bucket_offset: last_read_location.bucket_offset,
                        bucket_len: current_read_location.bucket_offset
                            + current_read_location.bucket_len
                            - last_read_location.bucket_offset,
                    },
                    from_star: current_star,
                    to_star: last_star,
                    from_exclamation_mark: current_exclamation_mark,
                    to_exclamation_mark: last_exclamation_mark,
                };

                let existing_edge = graph
                    .edges_between(reverse_from_node_index, reverse_to_node_index)
                    .next();
                if let Some(edge) = existing_edge {
                    let edge_data = graph.edge_data_mut(edge);
                    edge_data.add_edge_read_association(reverse_read_association);
                } /*else {
                      let edge_data = PlainWtdbg2EdgeData {
                          read_associations: vec![reverse_read_association],
                      };
                      graph.add_edge(
                          reverse_from_node_index,
                          reverse_to_node_index,
                          edge_data.into(),
                      );
                      warn!("Added non-existing edge N{} -> N{}", reverse_from_node_index.as_usize(), reverse_to_node_index.as_usize());
                  }*/
            }

            last_node = current_node;
        }
    }
    info!("Loaded {} edges", graph.edge_count());

    for edge_index in graph.edge_indices() {
        let edge_data = graph.edge_data(edge_index);
        if edge_data.edge_read_associations().is_empty() {
            warn!(
                "Edge has no read associations: N{} -> N{}",
                graph.edge_endpoints(edge_index).from_node.as_usize(),
                graph.edge_endpoints(edge_index).to_node.as_usize()
            );
        } else if edge_data.edge_read_associations().len() < 3 {
            let from_node = graph.edge_endpoints(edge_index).from_node.as_usize();
            let to_node = graph.edge_endpoints(edge_index).to_node.as_usize();
            warn!(
                "Edge has only {} read associations: N{}{} -> N{}{}",
                edge_data.edge_read_associations().len(),
                from_node / 2,
                if from_node % 2 == 1 { '+' } else { '-' },
                to_node / 2,
                if to_node % 2 == 1 { '+' } else { '-' }
            );
        }
    }

    Ok(graph)
}

/// Write a list of contigs in wtdbg's .ctg.lay format to a file.
pub fn write_contigs_to_wtdbg2_to_file<
    'ws,
    P1: AsRef<Path>,
    P2: AsRef<Path>,
    NodeData: Wtdbg2NodeData,
    EdgeData: Wtdbg2EdgeData,
    Graph: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph>,
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
    R: Read,
    W: Write,
    NodeData: Wtdbg2NodeData,
    EdgeData: Wtdbg2EdgeData,
    Graph: ImmutableGraphContainer<NodeData = NodeData, EdgeData = EdgeData>,
    Walk: 'ws + for<'w> EdgeWalk<'w, Graph>,
    WalkSource: 'ws + IntoIterator<Item = &'ws Walk>,
>(
    graph: &Graph,
    walks: WalkSource,
    raw_reads: bio::io::fasta::Reader<R>,
    output: &mut W,
) -> Result<()> {
    let mut read_map = HashMap::<_, VectorGenome>::new();
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
            genome.extend(record.seq().iter().copied());
        } else {
            read_map.insert(id.to_owned(), VectorGenome::from_iter(record.seq().iter()));
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
    walks.sort_by(|(l1, _), (l2, _)| l2.cmp(&l1));

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

        let first_edge_data = graph.edge_data(*walk.first().unwrap());
        let first_edge_length = first_edge_data.length();
        if first_edge_length < 4 {
            println!("First edge length: {}", first_edge_length);
            for read_association in first_edge_data.edge_read_associations() {
                println!("{:?}", read_association);
            }
        }
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
                    "S\t{}\t{}\t{}\t{}\t{}",
                    read_association.read_id,
                    if read_association.location.direction {
                        '+'
                    } else {
                        '-'
                    },
                    offset,
                    len,
                    VectorGenome::from_iter(
                        read_map.get(&read_association.read_id).unwrap()[offset..offset + len]
                            .iter()
                    )
                )?;
            }
        }
    }

    info!("Finished writing {} walks", printed_walks);
    info!("{} too short walks were dropped", dropped_walks);

    output.flush()?;
    Ok(())
}
