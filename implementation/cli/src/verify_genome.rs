use crate::CliOptions;
use clap::Clap;
use compact_genome::implementation::vec_sequence::AsciiVectorGenome;
use compact_genome::interface::sequence::GenomeSequence;

#[derive(Clap)]
pub struct VerifyGenomeCommand {
    #[clap(short, long, about = "The input file in fasta format")]
    pub input: String,
}

error_chain! {
    foreign_links {
        Io(std::io::Error);
    }

    errors {
        GenomeHasHole(invalid_characters: String) {
            description("the genome has a hole, i.e. its genome string contains characters other than ACGT")
            display("the genome has a hole")
        }

        GenomeHasNonUtf8Characters {
            description("the genome contains characters that are not valid UTF-8")
            display("the genome contains characters that are not valid UTF-8")
        }

        GenomeHasMultipleRecords {
            description("the genome consists of multiple fasta records")
            display("the genome consists of multiple fasta records")
        }

        GenomeHasNoRecords {
            description("the genome has no fasta records")
            display("the genome has no fasta records")
        }

        EmptyGenomeString {
            description("the genome string is empty")
            display("the genome string is empty")
        }
    }
}

pub(crate) fn verify_genome(
    _options: &CliOptions,
    subcommand: &VerifyGenomeCommand,
) -> crate::Result<()> {
    info!("Verifying that the genome has no holes...");

    info!("Reading genome from: {}", &subcommand.input);
    let records = bio::io::fasta::Reader::from_file(&subcommand.input)
        .map_err(|e| {
            error!("Error reading genome file");
            Error::from(e)
        })?
        .records();

    let mut records_found = 0;
    for record in records {
        records_found += 1;
        let genome: AsciiVectorGenome = match record {
            Ok(record) => record.seq().iter().copied().collect(),
            Err(err) => {
                error!("Error reading genome file");
                return Err(err.into());
            }
        };

        let invalid_characters = String::from_utf8(genome.get_invalid_characters());
        match invalid_characters {
            Ok(invalid_characters) => {
                if !invalid_characters.is_empty() {
                    error!(
                        "Genome contains a hole: invalid characters '{}'",
                        invalid_characters
                    );
                    return Err(Error::from(ErrorKind::GenomeHasHole(invalid_characters)).into());
                }
            }
            Err(_) => {
                error!("Genome contains a hole: characters that are not valid UTF-8");
                return Err(Error::from(ErrorKind::GenomeHasNonUtf8Characters).into());
            }
        }

        if genome.is_empty() {
            error!("Genome string is empty");
            return Err(Error::from(ErrorKind::EmptyGenomeString).into());
        }
    }

    if records_found == 0 {
        error!("Genome contains no fasta records");
        return Err(Error::from(ErrorKind::GenomeHasNoRecords).into());
    }

    info!("Found {} records", records_found);

    if records_found == 1 {
        info!("Genome has no holes");
    } else {
        info!("No record of the genome has a hole, but it consists of multiple records.")
    }

    Ok(())
}
