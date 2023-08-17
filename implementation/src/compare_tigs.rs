use crate::CliOptions;
use clap::Parser;
use compact_genome::implementation::bit_vec_sequence_store::BitVectorSequenceStore;
use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};
use compact_genome::interface::sequence::GenomeSequence;
use compact_genome::interface::sequence_store::SequenceStore;
use log::info;
use std::cmp::Ordering;
use std::collections::VecDeque;
use std::iter::{Enumerate, Peekable};
use std::mem;
use std::path::{Path, PathBuf};
use std::slice::Iter;
use suffix::SuffixTable;
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
    for (is_input2, index) in compare_sequence_files_by_equality(
        &sequences,
        |s1, s2, sequence_store| s1.cmp(s2, sequence_store),
        sequence_store,
    )? {
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

    for (file_index, sequences) in sequences.iter_mut().enumerate() {
        let file_index = file_index + 1;
        info!("Searching for duplicates within input{file_index}");

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

        info!("Removing duplicates from input{file_index}");
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

    info!("Building suffix array strings");
    let mut concatenated_strings = [String::new(), String::new()];
    let mut concatenated_string_offsets = [Vec::new(), Vec::new()];
    let string_terminator = '$';
    for index in 0..AlphabetType::CharacterType::ALPHABET_SIZE {
        let character = AlphabetType::CharacterType::from_index(index).unwrap();
        let character = AlphabetType::character_to_ascii(character);
        let character = char::from(character);
        assert_ne!(
            character, string_terminator,
            "Alphabet contains the string terminator {string_terminator:?}"
        );
    }

    for ((concatenated_string, sequences), concatenated_string_offsets) in concatenated_strings
        .iter_mut()
        .zip(sequences.iter())
        .zip(concatenated_string_offsets.iter_mut())
    {
        for sequence in sequences {
            concatenated_string.extend(sequence.iter(sequence_store).map(|alphabet_character| {
                char::from(AlphabetType::character_to_ascii(alphabet_character))
            }));
            concatenated_string.push(string_terminator);
            concatenated_string_offsets.push(concatenated_string.len());
        }
    }

    info!("Building suffix arrays");
    let suffix_arrays = [
        SuffixTable::new(&concatenated_strings[0]),
        SuffixTable::new(&concatenated_strings[1]),
    ];

    for (file_index, (sequences, (suffix_array, concatenated_string_offsets))) in sequences
        .iter_mut()
        .zip(suffix_arrays.iter().zip(concatenated_string_offsets.iter()))
        .enumerate()
    {
        let file_index = file_index + 1;
        info!("Searching for substrings within input{file_index}");

        let mut substrings = VecDeque::new();
        for (index1, sequence1) in sequences.iter().enumerate() {
            let sequence1_string: String = sequence1
                .iter(sequence_store)
                .map(|alphabet_character| {
                    char::from(AlphabetType::character_to_ascii(alphabet_character))
                })
                .collect();
            let sequence1_string_rc: String = sequence1
                .iter_reverse_complement(sequence_store)
                .map(|alphabet_character| {
                    char::from(AlphabetType::character_to_ascii(alphabet_character))
                })
                .collect();

            let mut positions: Vec<_> = suffix_array
                .positions(&sequence1_string)
                .iter()
                .copied()
                .map(|position| usize::try_from(position).unwrap())
                .collect();
            positions.extend(
                suffix_array
                    .positions(&sequence1_string_rc)
                    .iter()
                    .copied()
                    .map(|position| usize::try_from(position).unwrap()),
            );
            positions.sort_unstable();
            // println!("positions: {positions:?}");
            // println!("concatenated_string_offsets: {concatenated_string_offsets:?}");
            let mut indices2 = positions
                .into_iter()
                .map(|position| {
                    concatenated_string_offsets.partition_point(|&offset| offset <= position)
                })
                .fold(VecDeque::new(), |mut indices, position| {
                    if let Some(&back) = indices.back() {
                        if back != position {
                            indices.push_back(position);
                        }
                    } else {
                        indices.push_back(position);
                    }
                    indices
                });
            // println!("{indices2:?}");

            // Since we match against the whole suffix array we need to remove matches with ourselves.
            // Since we have sorted the sequences by length, the self-match is always the first.
            assert_eq!(
                indices2.front().copied(),
                Some(index1),
                "indices2: {indices2:?}"
            );
            indices2.pop_front();
            assert_ne!(
                indices2.front().copied(),
                Some(index1),
                "indices2: {indices2:?}"
            );

            if !indices2.is_empty() {
                substrings.push_back(index1);
            }

            for &index2 in &indices2 {
                assert!(
                    index2 < sequences.len(),
                    "indices2: {indices2:?}\nindex2: {index2}\nsequences.len(): {}",
                    sequences.len()
                );
                let sequence2 = &sequences[index2];
                info!(
                    "Input{}: sequence {} is substring of sequence {}",
                    file_index, sequence1.id, sequence2.id
                );
            }
        }

        info!("Removing substrings from input{file_index}");
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

    info!("Searching for sequences that are unique to one of the files");
    let unique = sequences;
    let mut unique1 = Vec::new();
    let mut unique2 = Vec::new();

    for (is_input2, index) in compare_sequence_files_by_equality(
        &unique,
        |s1, s2, sequence_store| match s1.len(sequence_store).cmp(&s2.len(sequence_store)) {
            Ordering::Equal => s1.cmp(s2, sequence_store),
            unequal => unequal,
        },
        sequence_store,
    )? {
        if !is_input2 {
            unique1.push(&unique[0][index]);
        } else {
            unique2.push(&unique[1][index]);
        }
    }
    /* println!(
        "unique1: {:?}",
        unique1
            .iter()
            .map(|sequence| &sequence.id)
            .collect::<Vec<_>>()
    );
    println!(
        "unique2: {:?}",
        unique2
            .iter()
            .map(|sequence| &sequence.id)
            .collect::<Vec<_>>()
    ); */
    let unique = [unique1, unique2];

    info!("Building suffix array strings for unique sequences");
    let mut concatenated_strings = [String::new(), String::new()];
    let mut concatenated_string_offsets = [Vec::new(), Vec::new()];

    for ((concatenated_string, unique), concatenated_string_offsets) in concatenated_strings
        .iter_mut()
        .zip(unique.iter())
        .zip(concatenated_string_offsets.iter_mut())
    {
        for sequence in unique {
            concatenated_string.extend(sequence.iter(sequence_store).map(|alphabet_character| {
                char::from(AlphabetType::character_to_ascii(alphabet_character))
            }));
            concatenated_string.push(string_terminator);
            concatenated_string_offsets.push(concatenated_string.len());
        }
    }

    info!("Building suffix arrays for unqiue sequences");
    let suffix_arrays = [
        SuffixTable::new(&concatenated_strings[0]),
        SuffixTable::new(&concatenated_strings[1]),
    ];

    for flipped in [false, true] {
        let (unique1, unique2, input1, input2) = if !flipped {
            (
                &unique[0],
                &unique[1],
                &subcommand.input1,
                &subcommand.input2,
            )
        } else {
            (
                &unique[1],
                &unique[0],
                &subcommand.input2,
                &subcommand.input1,
            )
        };

        info!(
            "Searching for substrings of input{} in input{}",
            if !flipped { 1 } else { 2 },
            if !flipped { 2 } else { 1 }
        );
        for sequence1 in unique1 {
            let sequence1_string: String = sequence1
                .iter(sequence_store)
                .map(|alphabet_character| {
                    char::from(AlphabetType::character_to_ascii(alphabet_character))
                })
                .collect();
            let sequence1_string_rc: String = sequence1
                .iter_reverse_complement(sequence_store)
                .map(|alphabet_character| {
                    char::from(AlphabetType::character_to_ascii(alphabet_character))
                })
                .collect();

            let (suffix_array, concatenated_string_offsets) = if !flipped {
                (&suffix_arrays[1], &concatenated_string_offsets[1])
            } else {
                (&suffix_arrays[0], &concatenated_string_offsets[0])
            };

            let mut positions: Vec<_> = suffix_array
                .positions(&sequence1_string)
                .iter()
                .copied()
                .map(|position| usize::try_from(position).unwrap())
                .collect();
            positions.extend(
                suffix_array
                    .positions(&sequence1_string_rc)
                    .iter()
                    .copied()
                    .map(|position| usize::try_from(position).unwrap()),
            );
            positions.sort_unstable();
            let indices2: Vec<_> = positions
                .into_iter()
                .map(|position| {
                    concatenated_string_offsets.partition_point(|&offset| offset <= position)
                })
                .fold(Vec::new(), |mut positions, position| {
                    if let Some(&last) = positions.last() {
                        if last != position {
                            positions.push(position);
                        }
                    } else {
                        positions.push(position);
                    }
                    positions
                });

            for index2 in indices2 {
                let sequence2 = &unique2[index2];

                assert!(
                    sequence1.len(sequence_store) < sequence2.len(sequence_store),
                    "sequence1 ({}) len: {}\n sequence2 ({}) len: {}",
                    sequence1.id,
                    sequence1.len(sequence_store),
                    sequence2.id,
                    sequence2.len(sequence_store)
                );

                info!(
                    "Input{} sequence {} is substring of input{} sequence {} ({:?}, {:?})",
                    if !flipped { 1 } else { 2 },
                    sequence1.id,
                    if !flipped { 2 } else { 1 },
                    sequence2.id,
                    input1,
                    input2
                );
            }
        }
    }

    Ok(())
}

struct CompareSortedSequencesByEquality<
    'a,
    AlphabetType: Alphabet,
    SequenceStoreType: SequenceStore<AlphabetType>,
    CompareFn: Fn(
        &DirectedSequence<AlphabetType, SequenceStoreType>,
        &DirectedSequence<AlphabetType, SequenceStoreType>,
        &SequenceStoreType,
    ) -> Ordering,
> {
    sequences1: Peekable<Enumerate<Iter<'a, DirectedSequence<AlphabetType, SequenceStoreType>>>>,
    sequences2: Peekable<Enumerate<Iter<'a, DirectedSequence<AlphabetType, SequenceStoreType>>>>,
    compare_fn: CompareFn,
    sequence_store: &'a SequenceStoreType,
}

impl<
        'a,
        AlphabetType: Alphabet,
        SequenceStoreType: SequenceStore<AlphabetType>,
        CompareFn: Fn(
            &DirectedSequence<AlphabetType, SequenceStoreType>,
            &DirectedSequence<AlphabetType, SequenceStoreType>,
            &SequenceStoreType,
        ) -> Ordering,
    > Iterator
    for CompareSortedSequencesByEquality<'a, AlphabetType, SequenceStoreType, CompareFn>
{
    type Item = (bool, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some((i1, s1)), Some((i2, s2))) =
            (self.sequences1.peek(), self.sequences2.peek())
        {
            match (self.compare_fn)(s1, s2, self.sequence_store) {
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
    CompareFn: 'a
        + Fn(
            &DirectedSequence<AlphabetType, SequenceStoreType>,
            &DirectedSequence<AlphabetType, SequenceStoreType>,
            &SequenceStoreType,
        ) -> Ordering,
>(
    sequences: &'a [Vec<DirectedSequence<AlphabetType, SequenceStoreType>>; 2],
    compare_fn: CompareFn,
    sequence_store: &'a SequenceStoreType,
) -> crate::Result<impl 'a + Iterator<Item = (bool, usize)>> {
    Ok(CompareSortedSequencesByEquality {
        sequences1: sequences[0].iter().enumerate().peekable(),
        sequences2: sequences[1].iter().enumerate().peekable(),
        compare_fn,
        sequence_store,
    })
}
