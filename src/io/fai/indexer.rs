use super::super::common;
use super::{ReadToFai, Record, Records};
use crate::errors::{Error, ErrorKind, Result};
use std::io::{BufRead, Seek};

/// Format represents the input format to be indexed
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Format {
    FASTA,
    FASTQ,
}

impl Format {
    /// Return the description prefix
    fn description_prefix(&self) -> u8 {
        match self {
            Format::FASTA => b'>',
            Format::FASTQ => b'@',
        }
    }

    /// Return the marker for the end of the sequence portion of the entry
    fn sequence_end_marker(&self) -> u8 {
        match self {
            Format::FASTA => self.description_prefix(),
            Format::FASTQ => b'+',
        }
    }
}

/// Indexer indexes input into Fai records
pub struct Indexer<R>
where
    R: std::io::Read + std::io::Seek,
{
    reader: std::io::BufReader<R>,
    format: Format,
    buffer: Vec<u8>,
    sequence_num_bytes: usize,
    eof: bool,
}

impl<R> Indexer<R>
where
    R: std::io::Read + std::io::Seek,
{
    /// Construct a new indexer
    pub fn new(reader: R, format: Format) -> Self {
        Self {
            reader: std::io::BufReader::new(reader),
            format,
            buffer: Vec::new(),
            sequence_num_bytes: 0,
            eof: false,
        }
    }

    /// Consume a reader by iterating over it
    pub fn iter(self) -> Records<Indexer<R>> {
        Records::new(self)
    }

    /// Read the first line of the input entry
    fn read_description(&mut self, record: &mut Record) -> Result<()> {
        if self.buffer.is_empty() {
            self.read_line()?;
        }
        record.name = get_name(&self.buffer, self.format)?;
        record.offset = self.reader.stream_position()?;
        self.buffer.clear();
        Ok(())
    }

    /// Read the entire sequence
    fn read_sequence(&mut self, record: &mut Record) -> Result<()> {
        self.sequence_num_bytes = 0;
        loop {
            if is_sequence_end(&self.buffer, self.format) || self.eof {
                return Ok(());
            }
            self.read_sequence_line(record)?;
        }
    }

    /// Read in a sequence line
    fn read_sequence_line(&mut self, record: &mut Record) -> Result<()> {
        self.buffer.clear();
        let num_bytes = self.read_line()?;
        if is_sequence_end(&self.buffer, self.format) {
            return Ok(());
        }
        if record.line_width == 0 {
            record.line_width = num_bytes;
            record.line_bases = common::count_bases(&self.buffer)?;
        } else if record.line_width < num_bytes {
            return Err(Error::new(ErrorKind::Input, "invalid record"));
        }
        self.sequence_num_bytes += num_bytes;
        record.length += common::count_bases(&self.buffer)?;
        Ok(())
    }

    /// Read in a line of data
    fn read_line(&mut self) -> Result<usize> {
        match common::read_line(&mut self.reader, &mut self.buffer) {
            Err(e) if e.kind == ErrorKind::Eof => {
                if self.eof {
                    Err(e)
                } else {
                    self.eof = true;
                    Ok(0)
                }
            }
            any => any,
        }
    }

    /// Read in the plus line
    fn read_plus(&mut self, record: &mut Record) -> Result<()> {
        if self.format == Format::FASTA {
            return Ok(());
        }
        record.qual_offset = Some(self.reader.stream_position()?);
        self.buffer.clear();
        Ok(())
    }

    /// Read the quality portion
    fn read_quality(&mut self) -> Result<()> {
        if self.format == Format::FASTA {
            return Ok(());
        }
        self.reader.consume(self.sequence_num_bytes);
        self.read_line()?;
        Ok(())
    }
}

impl<R> ReadToFai for Indexer<R>
where
    R: std::io::Read + std::io::Seek,
{
    /// Read a Fai record
    fn read(&mut self, record: &mut Record) -> Result<()> {
        self.read_description(record)?;
        self.read_sequence(record)?;
        self.read_plus(record)?;
        self.read_quality()?;
        Ok(())
    }
}

impl Indexer<std::fs::File> {
    /// Construct an indexer from path
    pub fn from_path<P: AsRef<std::path::Path>>(path: P, format: Format) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        Ok(Indexer::new(file, format))
    }
}

/// Check to see if the line is a description line
fn is_description(line: &[u8], format: Format) -> bool {
    line.starts_with(&[format.description_prefix()])
}

/// Retrieve name from the description line
fn get_name(description: &[u8], format: Format) -> Result<String> {
    if !is_description(description, format) {
        return Err(Error::new(ErrorKind::Input, "invalid input format"));
    }
    let description = String::from_utf8(description[1..].to_vec())?;
    Ok(common::parse_sequence_name(&description))
}

/// Check to see if the line has a sequence end marker
fn is_sequence_end(line: &[u8], format: Format) -> bool {
    line.starts_with(&[format.sequence_end_marker()])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_example() {
        // Test example from samtools documentation here:
        //   https://www.htslib.org/doc/faidx.html
        let input = br#">one
ATGCATGCATGCATGCATGCATGCATGCAT
GCATGCATGCATGCATGCATGCATGCATGC
ATGCAT
>two another chromosome
ATGCATGCATGCAT
GCATGCATGCATGC
"#;
        let input_line = std::io::Cursor::new(input);
        let mut indexer = Indexer::new(input_line, Format::FASTA);
        let mut record = Record::new();
        let expected = Record {
            name: "one".into(),
            length: 66,
            offset: 5,
            line_bases: 30,
            line_width: 31,
            qual_offset: None,
        };
        assert!(
            indexer.read(&mut record).is_ok(),
            "Should work for example in documentation",
        );
        assert_eq!(expected, record, "Should work for example in documentation");

        record.clear();
        let expected = Record {
            name: "two".into(),
            length: 28,
            offset: 98,
            line_bases: 14,
            line_width: 15,
            qual_offset: None,
        };
        assert!(
            indexer.read(&mut record).is_ok(),
            "Should read a second record",
        );
        assert_eq!(expected, record, "Should work for example in documentation");
    }

    #[test]
    fn test_fastq_example() {
        // Test example from samtools documentation here:
        //   https://www.htslib.org/doc/faidx.html
        let input = br#"@fastq1
ATGCATGCATGCATGCATGCATGCATGCAT
GCATGCATGCATGCATGCATGCATGCATGC
ATGCAT
+
FFFA@@FFFFFFFFFFHHB:::@BFFFFGG
HIHIIIIIIIIIIIIIIIIIIIIIIIFFFF
8011<<
@fastq2
ATGCATGCATGCAT
GCATGCATGCATGC
+
IIA94445EEII==
=>IIIIIIIIICCC"#;
        let input_line = std::io::Cursor::new(input);
        let mut reader = Indexer::new(input_line, Format::FASTQ);
        let mut record = Record::new();
        let expected = Record {
            name: "fastq1".into(),
            length: 66,
            offset: 8,
            line_bases: 30,
            line_width: 31,
            qual_offset: Some(79),
        };
        assert!(
            reader.read(&mut record).is_ok(),
            "Should work for example in documentation",
        );
        assert_eq!(expected, record, "Should work for example in documentation",);

        record.clear();
        let expected = Record {
            name: "fastq2".into(),
            length: 28,
            offset: 156,
            line_bases: 14,
            line_width: 15,
            qual_offset: Some(188),
        };
        assert!(
            reader.read(&mut record).is_ok(),
            "Should read a second record",
        );
        assert_eq!(expected, record, "Should work for example in documentation",);
    }

    #[test]
    fn test_is_description() {
        struct TestCase<'a> {
            name: &'a str,
            line: &'a [u8],
            format: Format,
            expected: bool,
        }
        let test_cases = [
            TestCase {
                name: "Should return true for a fasta description line",
                line: b">name",
                format: Format::FASTA,
                expected: true,
            },
            TestCase {
                name: "Should return true for a fastq description line",
                line: b"@name",
                format: Format::FASTQ,
                expected: true,
            },
            TestCase {
                name: "Should return false for a non description line",
                line: b"abcde",
                format: Format::FASTA,
                expected: false,
            },
            TestCase {
                name: "Should return false for an empty line",
                line: b"",
                format: Format::FASTQ,
                expected: false,
            },
        ];
        for test_case in test_cases {
            assert_eq!(
                test_case.expected,
                is_description(test_case.line, test_case.format),
                "{}",
                test_case.name,
            );
        }
    }

    #[test]
    fn test_get_name() {
        struct TestCase<'a> {
            name: &'a str,
            description: &'a [u8],
            format: Format,
            expect_error: bool,
            expected: String,
        }
        let test_cases = [
            TestCase {
                name: "Should return an error on an invalid description",
                description: b"abcdefg",
                format: Format::FASTQ,
                expect_error: true,
                expected: String::new(),
            },
            TestCase {
                name: "Should parse a name from the description for fasta formats",
                description: b">name bla",
                format: Format::FASTA,
                expect_error: false,
                expected: "name".to_string(),
            },
            TestCase {
                name: "Should parse a name from the description for fastq formats",
                description: b"@rita best pup",
                format: Format::FASTQ,
                expect_error: false,
                expected: "rita".to_string(),
            },
            TestCase {
                name: "Should skip spaces after the prefix before name",
                description: b"@     rita woof",
                format: Format::FASTQ,
                expect_error: false,
                expected: "rita".to_string(),
            },
        ];
        for test_case in test_cases {
            let actual = get_name(test_case.description, test_case.format);
            if test_case.expect_error {
                assert!(actual.is_err(), "{}", test_case.name);
            } else {
                assert_eq!(Ok(test_case.expected), actual, "{}", test_case.name);
            }
        }
    }

    #[test]
    fn test_is_sequence_end() {
        struct TestCase<'a> {
            name: &'a str,
            line: &'a [u8],
            format: Format,
            expected: bool,
        }
        let test_cases = [
            TestCase {
                name: "Should return true for a sequence end line for fasta",
                line: b">name",
                format: Format::FASTA,
                expected: true,
            },
            TestCase {
                name: "Should return true for a sequence end line for fastq",
                line: b"+abcde",
                format: Format::FASTQ,
                expected: true,
            },
            TestCase {
                name: "Should return false for a non sequence end line",
                line: b"AAGGCTT",
                format: Format::FASTA,
                expected: false,
            },
        ];
        for test_case in test_cases {
            assert_eq!(
                test_case.expected,
                is_sequence_end(test_case.line, test_case.format),
                "{}",
                test_case.name
            );
        }
    }
}
