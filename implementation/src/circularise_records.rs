use crate::CliOptions;
use bio::io::fasta::Record;
use clap::Parser;

#[derive(Parser)]
pub struct CirculariseGenomeCommand {
    #[clap(short, long, help = "The input file in fasta format")]
    pub input: String,

    #[clap(
        short,
        long,
        help = "The length of the prefix to append to the genome",
        default_value = "1000"
    )]
    pub circularisation_size: usize,

    #[clap(
        short,
        long,
        help = "The output file, to which the circularised records should be written as fasta"
    )]
    pub output: String,
}

pub(crate) fn circularise_records(
    _options: &CliOptions,
    subcommand: &CirculariseGenomeCommand,
) -> crate::Result<()> {
    info!(
        "Circularising each record of the genome with a circularisation size of {}",
        subcommand.circularisation_size
    );

    info!("Reading genome from: {}", &subcommand.input);
    let records = bio::io::fasta::Reader::from_file(&subcommand.input)
        .map_err(|e| {
            error!("Error reading genome file");
            e
        })?
        .records();
    info!("Creating/overwriting output file: {}", &subcommand.output);
    let mut writer = bio::io::fasta::Writer::to_file(&subcommand.output)?;
    let circularisation_size = subcommand.circularisation_size;

    let mut records_found = 0;
    for record in records {
        records_found += 1;
        let record = match record {
            Ok(record) => record,
            Err(err) => {
                error!("Error reading genome file");
                return Err(err.into());
            }
        };

        if record.seq().len() < circularisation_size {
            error!(
                "Record is shorter than circularisation size: {} < {}",
                record.seq().len(),
                circularisation_size
            );
            bail!("record is shorter than circularisation size");
        }

        let mut vec = Vec::from(record.seq());
        vec.extend_from_slice(&record.seq()[..circularisation_size]);
        let record = Record::with_attrs(record.id(), record.desc(), vec.as_slice());
        writer.write_record(&record)?;
    }

    if records_found == 0 {
        warn!("Genome contains no fasta records");
    } else {
        info!("Circularised {} records", records_found);
    }

    Ok(())
}
