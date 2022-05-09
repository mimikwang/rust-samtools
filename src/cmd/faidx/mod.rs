use crate::errors::{Error, ErrorKind, Result};
use crate::io::fai::{Fai, ReadToFai, Writer};
use crate::io::fasta::reader::fai::Reader as FastaReader;
use crate::io::fastq::reader::fai::Reader as FastqReader;

pub const SUBCOMMAND: &str = "faidx";
const FILE_ARG: &str = "file";
const FASTQ_FLAG: &str = "fastq";
const FASTQ_FLAG_SHORT: char = 'f';
const SUFFIX: &str = ".fai";

/// faidx subcommand
pub fn command() -> clap::Command<'static> {
    clap::Command::new(SUBCOMMAND)
        .arg(clap::Arg::new(FILE_ARG).required(true))
        .arg(
            clap::Arg::new(FASTQ_FLAG)
                .long(FASTQ_FLAG)
                .short(FASTQ_FLAG_SHORT)
                .takes_value(false),
        )
}

/// Run faidx workflow
pub fn run(matches: &clap::ArgMatches) -> Result<()> {
    let path = matches
        .value_of(FILE_ARG)
        .ok_or_else(|| Error::new(ErrorKind::User, "file argument required"))?;
    let mut reader = build_reader(path, matches)?;
    let mut writer = Writer::new(std::fs::File::create(output_name(path))?);
    loop {
        let mut record = Fai::new();
        let result = match reader.read(&mut record) {
            Ok(()) => Ok(record),
            Err(err) if err.kind == ErrorKind::Eof => {
                return Ok(());
            }
            Err(err) => Err(err),
        }?;
        writer.write(&result)?;
    }
}

/// Reader factory
fn build_reader(path: &str, matches: &clap::ArgMatches) -> Result<Box<dyn ReadToFai>> {
    if matches.is_present(FASTQ_FLAG) {
        return Ok(Box::new(FastqReader::from_path(path)?));
    }
    Ok(Box::new(FastaReader::from_path(path)?))
}

/// Output name for index file
fn output_name(file: &str) -> String {
    format!("{}{}", file, SUFFIX)
}
