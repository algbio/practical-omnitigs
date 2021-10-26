use crate::CliOptions;
use clap::Parser;
use std::io::Write;

#[derive(Parser)]
pub struct FilterCommand {
    #[clap(short, long, about = "The input file in fasta format")]
    pub input: String,

    #[clap(
        short,
        long,
        about = "The id of the record that should be retained, all other records are deleted"
    )]
    pub retain: Option<String>,

    #[clap(
        short,
        long,
        about = "The output file to which the remaining records should be written in fasta format"
    )]
    pub output: String,

    #[clap(
        short,
        long,
        about = "If given, extract the name of the genome into the given file"
    )]
    pub extract_name: Option<String>,
}

pub(crate) fn filter_records(
    _options: &CliOptions,
    subcommand: &FilterCommand,
) -> crate::Result<()> {
    if let Some(retain) = &subcommand.retain {
        info!("Retaining only the record with id '{}'", retain);
    } else {
        info!("Filtering nothing");
    }

    info!("Reading genome from: {}", &subcommand.input);
    let records = bio::io::fasta::Reader::from_file(&subcommand.input)
        .map_err(|e| {
            error!("Error reading genome file");
            e
        })?
        .records();
    info!("Creating/overwriting output file: {}", &subcommand.output);
    let mut writer = bio::io::fasta::Writer::to_file(&subcommand.output)?;
    let mut name_writer = if let Some(extract_name) = &subcommand.extract_name {
        info!("Creating/overwriting name file: {}", extract_name);
        Some(std::io::BufWriter::new(std::fs::File::create(
            extract_name,
        )?))
    } else {
        None
    };

    let mut records_read = 0;
    let mut records_written = 0;
    for record in records {
        records_read += 1;
        let record = match record {
            Ok(record) => record,
            Err(err) => {
                error!("Error reading genome file");
                return Err(err.into());
            }
        };

        let write_record = if let Some(retain) = &subcommand.retain {
            record.id() == retain
        } else {
            true
        };

        if write_record {
            writer.write_record(&record)?;
            records_written += 1;

            if let Some(name_writer) = &mut name_writer {
                writeln!(
                    name_writer,
                    "{} {}",
                    record.id(),
                    record.desc().unwrap_or("")
                )?;
            }
        }
    }

    if records_read == 0 {
        warn!("Genome contains no fasta records");
    } else if records_written == 0 {
        error!("Filtered out all records");
    } else if records_written == records_read {
        info!("Filtered out no records, {} are present", records_written);
    } else {
        info!(
            "Filtered out {} records, now {} remain",
            records_read - records_written,
            records_written
        );
    }

    Ok(())
}
