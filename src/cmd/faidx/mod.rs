use crate::errors::Result;
use crate::io::fai::Writer;
use crate::io::fasta::reader::fai::Reader as FastaReader;

pub const SUBCOMMAND: &str = "faidx";
const FILE_ARG: &str = "file";
const SUFFIX: &str = ".fai";

/// faidx subcommand
pub fn command() -> clap::Command<'static> {
    clap::Command::new(SUBCOMMAND).arg(clap::Arg::new(FILE_ARG).required(true))
}

/// Run faidx workflow
pub fn run(matches: &clap::ArgMatches) -> Result<()> {
    if let Some(file) = matches.value_of(FILE_ARG) {
        let reader = FastaReader::from_path(file)?;
        let output_file = std::fs::File::create(output_name(file))?;
        let mut writer = Writer::from_writer(output_file);
        for result in reader.iter() {
            let record = result?;
            writer.write(&record)?;
        }
    }
    Ok(())
}

/// Output name for index file
fn output_name(file: &str) -> String {
    format!("{}{}", file, SUFFIX)
}
