use crate::CliOptions;
use clap::Parser;
use compact_genome::implementation::bit_vec_sequence_store::BitVectorSequenceStore;
use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use compact_genome::interface::alphabet::Alphabet;
use compact_genome::interface::sequence::GenomeSequence;
use compact_genome::interface::sequence_store::SequenceStore;
use log::info;
use std::cmp::Ordering;
use std::collections::VecDeque;
use std::iter::{Enumerate, Peekable};
use std::mem;
use std::path::{Path, PathBuf};
use std::slice::Iter;
use traitsequence::interface::Sequence;

#[derive(Parser)]
pub struct CompareTigsCommand {
    #[clap(long)]
    /// Check if some sequences are substrings of others.
    pub check_substrings: bool,

    /// The first input fasta file.
    pub input1: PathBuf,

    /// The second input fasta file.
    pub input2: PathBuf,
}

struct DirectedSequence<AlphabetType: Alphabet, SequenceStoreType: SequenceStore<AlphabetType>> {
    pub id: String,
    pub handle: SequenceStoreType::Handle,
    pub forwards: bool,
}

impl<AlphabetType: Alphabet, SequenceStoreType: SequenceStore<AlphabetType>> Clone
    for DirectedSequence<AlphabetType, SequenceStoreType>
where
    SequenceStoreType::Handle: Clone,
{
    fn clone(&self) -> Self {
        Self {
            id: self.id.clone(),
            handle: self.handle.clone(),
            forwards: self.forwards,
        }
    }
}

impl<AlphabetType: Alphabet, SequenceStoreType: SequenceStore<AlphabetType>>
    DirectedSequence<AlphabetType, SequenceStoreType>
{
    pub fn canonicalise(&mut self, sequence_store: &SequenceStoreType) {
        if let Some((cf, cr)) = sequence_store
            .get(&self.handle)
            .iter()
            .cloned()
            .zip(sequence_store.get(&self.handle).reverse_complement_iter())
            .find(|(cf, cr)| cf != cr)
        {
            match cf.cmp(&cr) {
                Ordering::Less => self.forwards = true,
                Ordering::Equal => { /* ignore */ }
                Ordering::Greater => self.forwards = false,
            }
        }
    }

    pub fn iter<'this, 'sequence_store, 'result>(
        &'this self,
        sequence_store: &'sequence_store SequenceStoreType,
    ) -> Box<dyn 'result + DoubleEndedIterator<Item = AlphabetType::CharacterType>>
    where
        'this: 'result,
        'sequence_store: 'result,
    {
        if self.forwards {
            Box::new(sequence_store.get(&self.handle).iter().cloned())
        } else {
            Box::new(sequence_store.get(&self.handle).reverse_complement_iter())
        }
    }

    pub fn iter_reverse_complement<'this, 'sequence_store, 'result>(
        &'this self,
        sequence_store: &'sequence_store SequenceStoreType,
    ) -> Box<dyn 'result + DoubleEndedIterator<Item = AlphabetType::CharacterType>>
    where
        'this: 'result,
        'sequence_store: 'result,
    {
        if self.forwards {
            Box::new(sequence_store.get(&self.handle).reverse_complement_iter())
        } else {
            Box::new(sequence_store.get(&self.handle).iter().cloned())
        }
    }

    pub fn cmp(&self, other: &Self, sequence_store: &SequenceStoreType) -> Ordering {
        if let Some((c1, c2)) = self
            .iter(sequence_store)
            .zip(other.iter(sequence_store))
            .find(|(c1, c2)| c1 != c2)
        {
            c1.cmp(&c2)
        } else {
            let len1 = sequence_store.get(&self.handle).len();
            let len2 = sequence_store.get(&other.handle).len();
            len1.cmp(&len2)
        }
    }

    pub fn eq(&self, other: &Self, sequence_store: &SequenceStoreType) -> bool {
        let len1 = sequence_store.get(&self.handle).len();
        let len2 = sequence_store.get(&other.handle).len();
        if len1 == len2 {
            // comparing in reverse because we compare neighboring sorted strings
            self.iter(sequence_store)
                .rev()
                .zip(other.iter(sequence_store).rev())
                .all(|(c1, c2)| c1 == c2)
        } else {
            false
        }
    }

    pub fn len(&self, sequence_store: &SequenceStoreType) -> usize {
        sequence_store.get(&self.handle).len()
    }
}

fn read_fasta_sequences<AlphabetType: Alphabet, SequenceStoreType: SequenceStore<AlphabetType>>(
    file: impl AsRef<Path>,
    sequence_store: &mut SequenceStoreType,
) -> crate::Result<Vec<DirectedSequence<AlphabetType, SequenceStoreType>>> {
    let file = file.as_ref();
    info!("Reading fasta file {file:?}");
    let reader = bio::io::fasta::Reader::from_file(file)?;

    let mut result = Vec::new();
    for record in reader.records() {
        let record = record?;
        result.push(DirectedSequence {
            id: record.id().to_owned(),
            handle: sequence_store.add_from_slice_u8(record.seq())?,
            forwards: true,
        });
    }

    Ok(result)
}

pub(crate) fn compare_tigs(
    _options: &CliOptions,
    subcommand: &CompareTigsCommand,
) -> crate::Result<()> {
    let mut sequence_store = BitVectorSequenceStore::<DnaAlphabet>::new();

    info!("Loading sequences");
    let mut sequences = [
        read_fasta_sequences(&subcommand.input1, &mut sequence_store)?,
        read_fasta_sequences(&subcommand.input2, &mut sequence_store)?,
    ];

    info!("Canonicalising sequences");
    for sequences in &mut sequences {
        for sequence in sequences.iter_mut() {
            sequence.canonicalise(&sequence_store);
        }
    }

    if subcommand.check_substrings {
        compare_sequences_by_substrings(subcommand, sequences, &sequence_store)?;
    } else {
        compare_sequences_by_equality(subcommand, sequences, &sequence_store)?;
    }

    Ok(())
}

fn compare_sequences_by_equality<
    AlphabetType: Alphabet,
    SequenceStoreType: SequenceStore<AlphabetType>,
>(
    subcommand: &CompareTigsCommand,
    mut sequences: [Vec<DirectedSequence<AlphabetType, SequenceStoreType>>; 2],
    sequence_store: &SequenceStoreType,
) -> crate::Result<()> {
    info!("Sorting sequences alphabetically");
    for sequences in &mut sequences {
        sequences.sort_unstable_by(|s1, s2| s1.cmp(s2, sequence_store));
    }

    info!("Removing duplicates within files");
    for (file_index, sequences) in sequences.iter_mut().enumerate() {
        let file_index = file_index + 1;

        let mut duplicates = VecDeque::new();
        for (index, sequence_pair) in sequences.windows(2).enumerate() {
            if sequence_pair[0].eq(&sequence_pair[1], sequence_store) {
                duplicates.push_back(index + 1);
                info!(
                    "Input{}: sequence {} is equal to sequence {}",
                    file_index, sequence_pair[0].id, sequence_pair[1].id
                );
            }
        }

        let old_sequences = mem::take(sequences);
        sequences.extend(
            old_sequences
                .into_iter()
                .enumerate()
                .filter_map(|(index, sequence)| {
                    if let Some(first_duplicate) = duplicates.front().copied() {
                        if index == first_duplicate {
                            duplicates.pop_front();
                            None
                        } else {
                            Some(sequence)
                        }
                    } else {
                        Some(sequence)
                    }
                }),
        );
    }

    info!("Comparing sequence files");
    let sequences = sequences;
    for (is_input2, index) in compare_sequence_files_by_equality(&sequences, sequence_store)? {
        if !is_input2 {
            info!(
                "Sequence {} is unique to input1 ({:?})",
                sequences[0][index].id, subcommand.input1
            );
        } else {
            info!(
                "Sequence {} is unique to input2 ({:?})",
                sequences[1][index].id, subcommand.input2
            );
        }
    }

    Ok(())
}

fn compare_sequences_by_substrings<
    AlphabetType: Alphabet,
    SequenceStoreType: SequenceStore<AlphabetType>,
>(
    subcommand: &CompareTigsCommand,
    mut sequences: [Vec<DirectedSequence<AlphabetType, SequenceStoreType>>; 2],
    sequence_store: &SequenceStoreType,
) -> crate::Result<()> {
    info!("Sorting sequences by length, and sequences of the same length by alphabet");
    for sequences in &mut sequences {
        sequences.sort_unstable_by(|s1, s2| {
            match s1.len(sequence_store).cmp(&s2.len(sequence_store)) {
                Ordering::Equal => s1.cmp(s2, sequence_store),
                unequal => unequal,
            }
        });
    }

    info!("Removing substrings within files");
    for (file_index, sequences) in sequences.iter_mut().enumerate() {
        let file_index = file_index + 1;

        let mut substrings = VecDeque::new();
        for (index1, sequence1) in sequences.iter().enumerate() {
            let sequence1_vec: Vec<_> = sequence1.iter(sequence_store).collect();
            let sequence1_vec_rc: Vec<_> =
                sequence1.iter_reverse_complement(sequence_store).collect();

            'inner_loop: for sequence2 in sequences.iter().skip(index1 + 1) {
                if sequence1.len(sequence_store) == sequence2.len(sequence_store) {
                    if sequence1.eq(sequence2, sequence_store) {
                        substrings.push_back(index1);
                        info!(
                            "Input{}: sequence {} is equal to sequence {}",
                            file_index, sequence1.id, sequence2.id
                        );
                        break 'inner_loop;
                    }
                } else {
                    let sequence2_vec: Vec<_> = sequence2.iter(sequence_store).collect();

                    if sequence1_vec.is_proper_subsequence_of(&sequence2_vec)
                        || sequence1_vec_rc.is_proper_subsequence_of(&sequence2_vec)
                    {
                        substrings.push_back(index1);
                        info!(
                            "Input{}: sequence {} is substring of sequence {}",
                            file_index, sequence1.id, sequence2.id
                        );
                        break 'inner_loop;
                    }
                }
            }
        }

        let old_sequences = mem::take(sequences);
        sequences.extend(
            old_sequences
                .into_iter()
                .enumerate()
                .filter_map(|(index, sequence)| {
                    if let Some(first_substring) = substrings.front().copied() {
                        if index == first_substring {
                            substrings.pop_front();
                            None
                        } else {
                            Some(sequence)
                        }
                    } else {
                        Some(sequence)
                    }
                }),
        );
    }

    info!("Comparing sequence files");
    let sequences = sequences;
    let mut unique1 = Vec::new();
    let mut unique2 = Vec::new();

    for (is_input2, index) in compare_sequence_files_by_equality(&sequences, sequence_store)? {
        if !is_input2 {
            unique1.push(&sequences[0][index]);
        } else {
            unique2.push(&sequences[1][index]);
        }
    }

    for flipped in [false, true] {
        let (unique1, unique2, input1, input2) = if !flipped {
            (&unique1, &unique2, &subcommand.input1, &subcommand.input2)
        } else {
            (&unique2, &unique1, &subcommand.input2, &subcommand.input1)
        };

        for sequence1 in unique1 {
            let sequence1_vec: Vec<_> = sequence1.iter(sequence_store).collect();
            let sequence1_vec_rc: Vec<_> =
                sequence1.iter_reverse_complement(sequence_store).collect();

            let unique2_offset = unique2.partition_point(|sequence| {
                sequence1.len(sequence_store) <= sequence.len(sequence_store)
            });

            'inner_loop: for sequence2 in unique2.iter().skip(unique2_offset) {
                if sequence1.len(sequence_store) == sequence2.len(sequence_store) {
                    if sequence1.eq(sequence2, sequence_store) {
                        unreachable!("Equal sequences have been filtered before");
                    }
                } else {
                    let sequence2_vec: Vec<_> = sequence2.iter(sequence_store).collect();

                    if sequence1_vec.is_proper_subsequence_of(&sequence2_vec)
                        || sequence1_vec_rc.is_proper_subsequence_of(&sequence2_vec)
                    {
                        info!(
                            "Sequence {} in input{} is substring of sequence {} in input{} ({:?}, {:?})",
                            sequence1.id, if !flipped {1} else {2}, sequence2.id, if !flipped {2} else {1}, input1, input2
                        );
                        break 'inner_loop;
                    }
                }
            }
        }
    }

    Ok(())
}

struct CompareSortedSequencesByEquality<
    'a,
    AlphabetType: Alphabet,
    SequenceStoreType: SequenceStore<AlphabetType>,
> {
    sequences1: Peekable<Enumerate<Iter<'a, DirectedSequence<AlphabetType, SequenceStoreType>>>>,
    sequences2: Peekable<Enumerate<Iter<'a, DirectedSequence<AlphabetType, SequenceStoreType>>>>,
    sequence_store: &'a SequenceStoreType,
}

impl<'a, AlphabetType: Alphabet, SequenceStoreType: SequenceStore<AlphabetType>> Iterator
    for CompareSortedSequencesByEquality<'a, AlphabetType, SequenceStoreType>
{
    type Item = (bool, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some((i1, s1)), Some((i2, s2))) =
            (self.sequences1.peek(), self.sequences2.peek())
        {
            match s1.cmp(s2, self.sequence_store) {
                Ordering::Less => {
                    let i1 = *i1;
                    self.sequences1.next().unwrap();
                    return Some((false, i1));
                }
                Ordering::Equal => {
                    self.sequences1.next().unwrap();
                    self.sequences2.next().unwrap();
                }
                Ordering::Greater => {
                    let i2 = *i2;
                    self.sequences2.next().unwrap();
                    return Some((true, i2));
                }
            }
        }

        if let Some((i1, _)) = self.sequences1.next() {
            return Some((false, i1));
        }

        if let Some((i2, _)) = self.sequences2.next() {
            return Some((true, i2));
        }

        None
    }
}

fn compare_sequence_files_by_equality<
    'a,
    AlphabetType: Alphabet,
    SequenceStoreType: SequenceStore<AlphabetType>,
>(
    sequences: &'a [Vec<DirectedSequence<AlphabetType, SequenceStoreType>>; 2],
    sequence_store: &'a SequenceStoreType,
) -> crate::Result<impl 'a + Iterator<Item = (bool, usize)>> {
    Ok(CompareSortedSequencesByEquality {
        sequences1: sequences[0].iter().enumerate().peekable(),
        sequences2: sequences[1].iter().enumerate().peekable(),
        sequence_store,
    })
}
