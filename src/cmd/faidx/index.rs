use super::{get_file, output_name, FASTQ_FLAG};
use crate::errors::{ErrorKind, Result};
use crate::io::fai::{self, ReadToFai};
use crate::io::fasta::reader::fai::Reader as FastaReader;
use crate::io::fastq::reader::fai::Reader as FastqReader;

/// Run the indexing workflow
pub fn run(matches: &clap::ArgMatches) -> Result<()> {
    let file = get_file(matches)?;
    let fastq = matches.is_present(FASTQ_FLAG);
    let mut reader = build_reader(file, fastq)?;
    let mut writer = fai::Writer::new(std::fs::File::create(output_name(file))?);
    consume_reader(&mut reader, &mut writer)?;
    Ok(())
}

/// Build the appropriate Fai reader
///
/// If `fastq` is true, then return a FASTQ Fai record reader.  Otherwise, return a FASTA Fai
/// record.
///
fn build_reader(file: &str, fastq: bool) -> Result<Box<dyn ReadToFai>> {
    if fastq {
        Ok(Box::new(FastqReader::from_path(file)?))
    } else {
        Ok(Box::new(FastaReader::from_path(file)?))
    }
}

/// Consume a reader and write to output
///
/// Duplicate sequence names are ignored.
///
fn consume_reader<W>(reader: &mut Box<dyn ReadToFai>, writer: &mut fai::Writer<W>) -> Result<()>
where
    W: std::io::Write,
{
    let mut names = std::collections::HashSet::<String>::new();
    loop {
        let mut record = fai::Record::new();
        let result = match reader.read(&mut record) {
            Ok(()) => Ok(record),
            Err(err) if err.kind == ErrorKind::Eof => {
                return Ok(());
            }
            Err(err) => Err(err),
        }?;
        if names.contains(&result.name) {
            println!("Duplicate entry: {}, skipping", &result.name)
        } else {
            writer.write(&result)?;
            names.insert(result.name);
        }
    }
}
