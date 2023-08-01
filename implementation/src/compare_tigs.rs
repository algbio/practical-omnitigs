use crate::CliOptions;
use clap::Parser;
use compact_genome::implementation::bit_vec_sequence_store::BitVectorSequenceStore;
use compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use compact_genome::interface::alphabet::Alphabet;
use compact_genome::interface::sequence::GenomeSequence;
use compact_genome::interface::sequence_store::SequenceStore;
use log::info;
use std::cmp::Ordering;
use std::path::{Path, PathBuf};
use traitsequence::interface::Sequence;

#[derive(Parser)]
pub struct CompareTigsCommand {
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
    ) -> Box<dyn 'result + Iterator<Item = AlphabetType::CharacterType>>
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

    info!("Sorting sequences");
    for sequences in &mut sequences {
        sequences.sort_unstable_by(|s1, s2| s1.cmp(s2, &sequence_store));
    }

    info!("Comparing sequences");
    let mut sequences1 = sequences[0].iter().peekable();
    let mut sequences2 = sequences[1].iter().peekable();

    while let (Some(s1), Some(s2)) = (sequences1.peek(), sequences2.peek()) {
        match s1.cmp(s2, &sequence_store) {
            Ordering::Less => {
                info!(
                    "Sequence {} is unique to input1 ({:?})",
                    s1.id, subcommand.input1
                );
                sequences1.next().unwrap();
            }
            Ordering::Equal => {
                sequences1.next().unwrap();
                sequences2.next().unwrap();
            }
            Ordering::Greater => {
                info!(
                    "Sequence {} is unique to input2 ({:?})",
                    s2.id, subcommand.input2
                );
                sequences2.next().unwrap();
            }
        }
    }

    for s1 in sequences1 {
        info!(
            "Sequence {} is unique to input1 ({:?})",
            s1.id, subcommand.input1
        );
    }

    for s2 in sequences2 {
        info!(
            "Sequence {} is unique to input2 ({:?})",
            s2.id, subcommand.input2
        );
    }

    Ok(())
}
