use crate::errors::{Error, ErrorKind, Result};

mod faidx;

/// Run the command line
pub fn run() -> Result<()> {
    let matches = clap::Command::new("rust-samtools")
        .author("Mimi Wang, mimikwang@gmail.com")
        .version("0.1.0")
        .about("Rust implementation of samtools")
        .subcommand(faidx::command())
        .subcommand_required(true)
        .get_matches();

    match matches.subcommand() {
        Some((faidx::SUBCOMMAND, matches)) => faidx::run(matches),
        Some((subcommand, _)) => Err(Error::new(
            ErrorKind::User,
            &format!("unrecognized command {}", subcommand),
        )),
        None => Ok(()),
    }
}
