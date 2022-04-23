mod reader;

use crate::errors::{Error, ErrorKind, Result};
use serde::{Deserialize, Serialize};

pub use reader::{IterFai, Reader};

const DELIMITER: u8 = b'\t';

/// An fai entry - description at the following:
///     https://www.htslib.org/doc/faidx.html
#[derive(Debug, Deserialize, PartialEq, Serialize, Default)]
pub struct Fai {
    /// Name of this reference sequence
    pub name: String,
    /// Total length of this reference sequence, in bases
    pub length: usize,
    /// Offset in the FASTA/FASTQ file of this sequence's first base
    pub offset: usize,
    /// The number of bases on each line
    pub line_bases: usize,
    /// The number of bytes in each line, including the newline
    pub line_width: usize,
    /// Offset of sequence's first quality within the FASTQ file
    pub qual_offset: Option<usize>,
}

impl Fai {
    /// Construct a default fai
    pub fn new() -> Self {
        Self::default()
    }
}

impl std::convert::TryFrom<&mut csv::StringRecord> for Fai {
    type Error = Error;

    fn try_from(record: &mut csv::StringRecord) -> Result<Self> {
        if record.len() > 6 {
            return Err(Error::new(ErrorKind::Input, "invalid fai format"));
        }
        if record.len() == 5 {
            record.push_field("");
        }
        let fai: Fai = record.deserialize(None)?;
        Ok(fai)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fai_from_string_record() {
        struct TestCase<'a> {
            name: &'a str,
            record: csv::StringRecord,
            expect_error: bool,
            expected: Fai,
        }
        let test_cases = &mut [
            TestCase {
                name: "Should parse a valid string record with 6 entries [Issue #5]",
                record: csv::StringRecord::from(vec!["name", "1", "1", "10", "4", "3"]),
                expect_error: false,
                expected: Fai {
                    name: "name".into(),
                    length: 1,
                    offset: 1,
                    line_bases: 10,
                    line_width: 4,
                    qual_offset: Some(3),
                },
            },
            TestCase {
                name: "Should parse a valid string record with 5 entries [Issue #5]",
                record: csv::StringRecord::from(vec!["name", "1", "1", "10", "4"]),
                expect_error: false,
                expected: Fai {
                    name: "name".into(),
                    length: 1,
                    offset: 1,
                    line_bases: 10,
                    line_width: 4,
                    qual_offset: None,
                },
            },
            TestCase {
                name: "Should return an error with less than 5 entries [Issue #5]",
                record: csv::StringRecord::from(vec!["name", "1", "1", "10"]),
                expect_error: true,
                expected: Fai::new(),
            },
            TestCase {
                name: "Should return an error with more than 6 entries [Issue #5]",
                record: csv::StringRecord::from(vec!["name", "1", "1", "10", "4", "6", "6"]),
                expect_error: true,
                expected: Fai::new(),
            },
            TestCase {
                name: "Should return an error if entries are not the right type [Issue #5]",
                record: csv::StringRecord::from(vec!["name", "1", "1.2", "10", "4", "6"]),
                expect_error: true,
                expected: Fai::new(),
            },
        ];

        for test_case in test_cases.iter_mut() {
            let actual = Fai::try_from(&mut test_case.record);
            if test_case.expect_error {
                assert!(actual.is_err(), "{}", test_case.name);
            } else {
                assert!(actual.is_ok(), "{}", test_case.name);
                let actual = actual.unwrap();
                assert_eq!(test_case.expected, actual, "{}", test_case.name);
            }
        }
    }
}
