use crate::CliOptions;
use clap::Clap;

#[derive(Clap)]
pub struct FilterCommand {
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
}

pub(crate) fn filter_records(
    options: &CliOptions,
    subcommand: &FilterCommand,
) -> crate::Result<()> {
    if let Some(retain) = &subcommand.retain {
        info!("Retaining only the record with id '{}'", retain);
    } else {
        info!("Filtering nothing");
    }

    info!("Reading genome from: {}", &options.input);
    let records = bio::io::fasta::Reader::from_file(&options.input)
        .map_err(|e| {
            error!("Error reading genome file");
            e
        })?
        .records();
    info!("Creating/overwriting output file: {}", &subcommand.output);
    let mut writer = bio::io::fasta::Writer::to_file(&subcommand.output)?;

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

        if let Some(retain) = &subcommand.retain {
            if record.id() == retain {
                writer.write_record(&record)?;
                records_written += 1;
            }
        } else {
            writer.write_record(&record)?;
            records_written += 1;
        }
    }

    if records_read == 0 {
        warn!("Genome contains no fasta records");
    } else if records_written == 0 {
        error!("Filtered out all records");
    } else {
        info!(
            "Filtered out {} records, now {} remain",
            records_read - records_written,
            records_written
        );
    }

    Ok(())
}
