use super::*;
use crate::errors::{Error, ErrorKind, Result};
use crate::io::fai::{Fai, IterFai, ReadToFai};

/// Reader reads a FASTA file into Fai records
pub struct Reader<B>
where
    B: std::io::BufRead + std::io::Seek,
{
    reader: B,
    buffer: Vec<u8>,
    eof: bool,
}

impl<B> Reader<B>
where
    B: std::io::BufRead + std::io::Seek,
{
    /// Construct a new reader
    pub fn new(reader: B) -> Self {
        Self {
            reader,
            buffer: Vec::new(),
            eof: false,
        }
    }

    /// Consume the reader and return an `IterFai`
    pub fn iter(self) -> IterFai<Self> {
        IterFai::new(self)
    }

    /// Read the first line of the FASTA entry
    fn read_description(&mut self, record: &mut Fai) -> Result<()> {
        if self.buffer.is_empty() {
            self.read_line()?;
        };
        record.name = get_name(&self.buffer)?;
        record.offset = self.reader.stream_position()?;
        self.buffer.clear();
        Ok(())
    }

    /// Read the entire sequence
    fn read_sequence(&mut self, record: &mut Fai) -> Result<()> {
        loop {
            if is_description(&self.buffer) || self.eof {
                return Ok(());
            }
            self.read_sequence_line(record)?;
        }
    }

    /// Read in a sequence line
    fn read_sequence_line(&mut self, record: &mut Fai) -> Result<()> {
        self.buffer.clear();
        let num_bytes = self.read_line()?;
        if is_description(&self.buffer) {
            return Ok(());
        }
        if record.line_width == 0 {
            record.line_width = num_bytes;
            record.line_bases = common::count_bases(&self.buffer)?;
        } else if record.line_width < num_bytes {
            return Err(Error::new(ErrorKind::Input, "invalid fasta record"));
        }
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
}

impl<R> ReadToFai for Reader<R>
where
    R: std::io::BufRead + std::io::Seek,
{
    /// Read a Fai record
    fn read(&mut self, record: &mut Fai) -> Result<()> {
        self.read_description(record)?;
        self.read_sequence(record)?;
        Ok(())
    }
}

impl<R> Reader<std::io::BufReader<R>>
where
    R: std::io::Read + std::io::Seek,
{
    /// Construct from a read seeker
    pub fn from_readseeker(readseeker: R) -> Self {
        Reader::new(std::io::BufReader::new(readseeker))
    }
}

impl Reader<std::io::BufReader<std::fs::File>> {
    /// Construct a reader from path
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        Ok(Reader::from_readseeker(file))
    }
}

/// Retrieve name from the description line
fn get_name(description: &[u8]) -> Result<String> {
    if !is_description(description) {
        return Err(Error::new(ErrorKind::Input, "invalid fasta format"));
    }
    let description = String::from_utf8(description[1..].to_vec())?;
    Ok(common::parse_sequence_name(&description))
}

/// Check to see if the line is a description line
fn is_description(line: &[u8]) -> bool {
    line.starts_with(&[SEQ_START])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_example_from_samtools_doc() {
        // Test example from samtools documentation here:
        //  https://www.htslib.org/doc/faidx.html
        let input = br#">one
ATGCATGCATGCATGCATGCATGCATGCAT
GCATGCATGCATGCATGCATGCATGCATGC
ATGCAT
>two another chromosome
ATGCATGCATGCAT
GCATGCATGCATGC"#;
        let input_line = std::io::Cursor::new(input);
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        let expected = Fai {
            name: "one".into(),
            length: 66,
            offset: 5,
            line_bases: 30,
            line_width: 31,
            qual_offset: None,
        };
        assert!(
            reader.read(&mut record).is_ok(),
            "Should work for example in documentation",
        );
        assert_eq!(expected, record, "Should work for example in documentation",);

        record.clear();
        let expected = Fai {
            name: "two".into(),
            length: 28,
            offset: 98,
            line_bases: 14,
            line_width: 15,
            qual_offset: None,
        };
        assert!(
            reader.read(&mut record).is_ok(),
            "Should read a second record",
        );
        assert_eq!(expected, record, "Should work for example in documentation",);
    }

    #[test]
    fn test_reader_read() {
        let input_line = std::io::Cursor::new(b">abc aa\nAAAA\nAAA\n>abcdef\nAAAAA\n");
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        let expected = Fai {
            name: "abc".into(),
            offset: 8,
            line_width: 5,
            line_bases: 4,
            length: 7,
            qual_offset: None,
        };
        assert!(
            reader.read(&mut record).is_ok(),
            "Should read FASTA entry into Fai record",
        );
        assert_eq!(expected, record, "Should read in a record");

        record.clear();
        let expected = Fai {
            name: "abcdef".into(),
            offset: 25,
            line_width: 6,
            line_bases: 5,
            length: 5,
            qual_offset: None,
        };
        assert!(
            reader.read(&mut record).is_ok(),
            "Should read a second record",
        );
        assert_eq!(expected, record, "{}", "Should read in a second record",);

        record.clear();
        assert_eq!(
            ErrorKind::Eof,
            reader.read(&mut record).unwrap_err().kind,
            "Should return an Eof error",
        );
    }

    #[test]
    fn test_reader_read_description() {
        let input_line = std::io::Cursor::new(b">abcdef aaaa\nAAAA\nAAAA\n");
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        assert!(
            reader.read_description(&mut record).is_ok(),
            "Should read description into Fai record",
        );
        assert_eq!("abcdef", &record.name, "Should read in the record name",);
        assert_eq!(13, record.offset, "Should read in the sequence offset",);
    }

    #[test]
    fn test_reader_read_sequence() {
        let input_line = std::io::Cursor::new(b"AAAA\nAAAA\n>abcdef\nAAAA\n");
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        assert!(
            reader.read_sequence(&mut record).is_ok(),
            "Should read sequence into Fai record"
        );
        assert_eq!(
            4, record.line_bases,
            "Should read in the number of bases in a line",
        );
        assert_eq!(
            5, record.line_width,
            "Should read in the number of bytes in a line",
        );
        assert_eq!(8, record.length, "Should increment number of bases",);

        let input_line = std::io::Cursor::new(b"AAAA\nAAAAA\n>abcdef\nAAAA\n");
        let mut reader = Reader::new(input_line);
        record.clear();
        assert!(
            reader.read_sequence(&mut record).is_err(),
            "Should error out if the line widths are inconsistent"
        );
    }

    #[test]
    fn test_reader_read_sequence_line() {
        let input_line = std::io::Cursor::new(b"AAAA\r\nAAAA\r\n");
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        assert!(
            reader.read_sequence_line(&mut record).is_ok(),
            "Should read sequence line into Fai record"
        );
        assert_eq!(
            4, record.line_bases,
            "Should read in the number of bases in a line",
        );
        assert_eq!(
            6, record.line_width,
            "Should read in the number of bytes in a line",
        );
        assert_eq!(4, record.length, "Should increment number of bases",);
    }

    #[test]
    fn test_get_name() {
        struct TestCase<'a> {
            name: &'a str,
            description: &'a [u8],
            expect_error: bool,
            expected: String,
        }
        let test_cases = [
            TestCase {
                name: "Should return an error on an invalid description",
                description: b"abcdefg",
                expect_error: true,
                expected: String::new(),
            },
            TestCase {
                name: "Should parse a name from the description",
                description: b">name a b c d e f",
                expect_error: false,
                expected: "name".to_string(),
            },
            TestCase {
                name: "Should skip spaces after > before name",
                description: b">    name abcdef",
                expect_error: false,
                expected: "name".to_string(),
            },
        ];
        for test_case in test_cases {
            let actual = get_name(test_case.description);
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
    fn test_is_description() {
        struct TestCase<'a> {
            name: &'a str,
            line: &'a [u8],
            expected: bool,
        }
        let test_cases = [
            TestCase {
                name: "Should return true for a description line",
                line: b">abcdef",
                expected: true,
            },
            TestCase {
                name: "Should return false for a non description line",
                line: b"bla",
                expected: false,
            },
            TestCase {
                name: "Should return false for an empty line",
                line: b"",
                expected: false,
            },
        ];
        for test_case in test_cases {
            assert_eq!(
                test_case.expected,
                is_description(test_case.line),
                "{}",
                test_case.name
            );
        }
    }

    #[test]
    fn test_reader_iter() {
        let input_line = std::io::Cursor::new(b">abc aa\nAAAA\nAAA\n>abcdef\nAAAAA\n");
        let results: Vec<Result<Fai>> = Reader::new(input_line).iter().collect();
        assert_eq!(2, results.len(), "Should iterate through all records");
    }
}
