use super::{get_file, output_name, FASTQ_FLAG};
use crate::errors::Result;
use crate::io::fai;
use std::fs::File;

/// Run the indexing workflow
pub fn run(matches: &clap::ArgMatches) -> Result<()> {
    let file = get_file(matches)?;
    let fastq = matches.is_present(FASTQ_FLAG);
    let reader = build_reader(file, fastq)?;
    let mut writer = fai::Writer::new(File::create(output_name(file))?);
    consume_reader(reader, &mut writer)
}

/// Build the appropriate Fai reader
///
/// If `fastq` is true, then return a FASTQ Fai record reader.  Otherwise, return a FASTA Fai
/// record.
///
fn build_reader(file: &str, fastq: bool) -> Result<fai::Indexer<File>> {
    let format = if fastq {
        fai::IndexerFormat::FASTQ
    } else {
        fai::IndexerFormat::FASTA
    };
    fai::Indexer::from_path(file, format)
}

/// Consume a reader and write to output
///
/// Duplicate sequence names are ignored.
///
fn consume_reader<W>(reader: fai::Indexer<File>, writer: &mut fai::Writer<W>) -> Result<()>
where
    W: std::io::Write,
{
    let mut names = std::collections::HashSet::<String>::new();
    for result in reader.iter() {
        let record = result?;
        if names.contains(&record.name) {
            println!("duplicate entry: {}, skipping", &record.name);
            continue;
        }
        writer.write(&record)?;
        names.insert(record.name);
    }
    Ok(())
}
