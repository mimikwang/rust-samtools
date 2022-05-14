mod indexer;
mod reader;
mod writer;

use crate::errors::{Error, ErrorKind, Result};
use serde::{Deserialize, Serialize};

pub use indexer::{Format as IndexerFormat, Indexer};
pub use reader::Reader;
pub use writer::Writer;

const DELIMITER: u8 = b'\t';
const FASTA_WIDTH: usize = 5;
const FASTQ_WIDTH: usize = 6;

/// ReadToFai reads to a Fai record
pub trait ReadToFai {
    fn read(&mut self, record: &mut Record) -> Result<()>;
}

/// Fai record as defined in the [`documentation`]
///
/// [`documentation`]: https://www.htslib.org/doc/faidx.html
///
#[derive(Debug, Deserialize, PartialEq, Serialize, Default)]
pub struct Record {
    /// Name of this reference sequence
    pub name: String,
    /// Total length of this reference sequence, in bases
    pub length: usize,
    /// Offset in the FASTA/FASTQ file of this sequence's first base
    pub offset: u64,
    /// The number of bases on each line
    pub line_bases: usize,
    /// The number of bytes in each line, including the newline
    pub line_width: usize,
    /// Offset of sequence's first quality within the FASTQ file.  This is always `None` for
    /// FASTA fies.
    pub qual_offset: Option<u64>,
}

impl Record {
    /// Construct a default Fai record
    pub fn new() -> Self {
        Self::default()
    }

    /// Convert a Fai record to a csv string record
    pub fn to_string_record(&self) -> csv::StringRecord {
        let mut record = vec![
            self.name.clone(),
            self.length.to_string(),
            self.offset.to_string(),
            self.line_bases.to_string(),
            self.line_width.to_string(),
        ];
        if let Some(qual_offset) = self.qual_offset {
            record.push(qual_offset.to_string());
        }
        csv::StringRecord::from(record)
    }

    /// Clear the Fai record
    pub fn clear(&mut self) {
        self.name.clear();
        self.length = 0;
        self.offset = 0;
        self.line_bases = 0;
        self.line_width = 0;
        self.qual_offset = None;
    }
}

impl std::convert::TryFrom<&mut csv::StringRecord> for Record {
    type Error = Error;

    /// Convert a csv string record to a Fai record
    ///
    /// An error is returned if it doesn't conform to either the FASTA or FASTQ formats.
    ///
    fn try_from(record: &mut csv::StringRecord) -> Result<Self> {
        if record.len() > FASTQ_WIDTH {
            return Err(Error::new(ErrorKind::Input, "invalid fai format"));
        }
        if record.len() == FASTA_WIDTH {
            record.push_field("");
        }
        let record: Record = record.deserialize(None)?;
        Ok(record)
    }
}

/// Type for iterating over fai records
#[derive(Debug)]
pub struct Records<F>
where
    F: ReadToFai,
{
    reader: F,
}

impl<F> Records<F>
where
    F: ReadToFai,
{
    /// Construct a new Records given a reader
    pub fn new(reader: F) -> Self {
        Self { reader }
    }
}

impl<F> Iterator for Records<F>
where
    F: ReadToFai,
{
    type Item = Result<Record>;

    /// Implementation of the Iterator trait for IterFai
    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Ok(()) => Some(Ok(record)),
            Err(err) if err.kind == ErrorKind::Eof => None,
            Err(err) => Some(Err(err)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fai_try_from_string_record() {
        struct TestCase<'a> {
            name: &'a str,
            record: csv::StringRecord,
            expect_error: bool,
            expected: Record,
        }
        let test_cases = &mut [
            TestCase {
                name: "Should parse a valid string record with 6 entries",
                record: csv::StringRecord::from(vec!["name", "1", "1", "10", "4", "3"]),
                expect_error: false,
                expected: Record {
                    name: "name".into(),
                    length: 1,
                    offset: 1,
                    line_bases: 10,
                    line_width: 4,
                    qual_offset: Some(3),
                },
            },
            TestCase {
                name: "Should parse a valid string record with 5 entries",
                record: csv::StringRecord::from(vec!["name", "1", "1", "10", "4"]),
                expect_error: false,
                expected: Record {
                    name: "name".into(),
                    length: 1,
                    offset: 1,
                    line_bases: 10,
                    line_width: 4,
                    qual_offset: None,
                },
            },
            TestCase {
                name: "Should return an error with less than 5 entries",
                record: csv::StringRecord::from(vec!["name", "1", "1", "10"]),
                expect_error: true,
                expected: Record::new(),
            },
            TestCase {
                name: "Should return an error with more than 6 entries",
                record: csv::StringRecord::from(vec!["name", "1", "1", "10", "4", "6", "6"]),
                expect_error: true,
                expected: Record::new(),
            },
            TestCase {
                name: "Should return an error if entries are not the right type",
                record: csv::StringRecord::from(vec!["name", "1", "1.2", "10", "4", "6"]),
                expect_error: true,
                expected: Record::new(),
            },
        ];

        for test_case in test_cases.iter_mut() {
            let actual = Record::try_from(&mut test_case.record);
            if test_case.expect_error {
                assert!(actual.is_err(), "{}", test_case.name);
            } else {
                assert!(actual.is_ok(), "{}", test_case.name);
                let actual = actual.unwrap();
                assert_eq!(test_case.expected, actual, "{}", test_case.name);
            }
        }
    }

    #[test]
    fn test_fai_to_string_record() {
        struct TestCase<'a> {
            name: &'a str,
            record: Record,
            expected: csv::StringRecord,
        }
        let test_cases = &[
            TestCase {
                name: "Should convert a fasta fai record",
                record: Record {
                    name: "name".into(),
                    length: 1,
                    offset: 2,
                    line_bases: 3,
                    line_width: 4,
                    qual_offset: None,
                },
                expected: csv::StringRecord::from(vec!["name", "1", "2", "3", "4"]),
            },
            TestCase {
                name: "Should convert a fastq fai record",
                record: Record {
                    name: "name".into(),
                    length: 1,
                    offset: 2,
                    line_bases: 3,
                    line_width: 4,
                    qual_offset: Some(5),
                },
                expected: csv::StringRecord::from(vec!["name", "1", "2", "3", "4", "5"]),
            },
        ];

        for test_case in test_cases.iter() {
            let actual = test_case.record.to_string_record();
            assert_eq!(test_case.expected, actual, "{}", test_case.name);
        }
    }
}
