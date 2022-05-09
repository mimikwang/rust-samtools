//! Rust implementation of [`samtools`]
//!
//! [`samtools`]: https://www.htslib.org/
mod cmd;
mod errors;
mod io;

extern crate clap;
extern crate csv;
extern crate serde;

fn main() -> errors::Result<()> {
    cmd::run()
}
