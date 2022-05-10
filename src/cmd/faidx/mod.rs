use crate::errors::{Error, ErrorKind, Result};

mod index;

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
    index::run(matches)
}

/// Get file argument
fn get_file(matches: &clap::ArgMatches) -> Result<&str> {
    matches
        .value_of(FILE_ARG)
        .ok_or_else(|| Error::new(ErrorKind::User, "file argument required"))
}

/// Output name for index file
fn output_name(file: &str) -> String {
    format!("{}{}", file, SUFFIX)
}
