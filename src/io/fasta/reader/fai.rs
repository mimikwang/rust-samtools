use super::*;
use crate::errors::{Error, ErrorKind, Result};
use crate::io::fai::Fai;

/// Reader reads a fasta file into FAI entries
pub struct Reader<B>
where
    B: std::io::BufRead + std::io::Seek,
{
    reader: B,
    buffer: Vec<u8>,
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
        }
    }

    /// Read a FAI record
    pub fn read(&mut self, record: &mut Fai) -> Result<()> {
        self.read_description(record)?;
        self.read_sequence(record)?;
        Ok(())
    }

    /// Read the first line of the fasta entry
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
            if is_description(&self.buffer) {
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
            record.line_bases = count_bases(&self.buffer)?;
        } else if record.line_width < num_bytes {
            return Err(Error::new(ErrorKind::Input, "invalid fasta record"));
        }
        record.length += count_bases(&self.buffer)?;
        Ok(())
    }

    /// Read in a line of data
    fn read_line(&mut self) -> Result<usize> {
        let num_bytes = self.reader.read_until(NEWLINE, &mut self.buffer)?;
        if num_bytes == 0 {
            return Err(Error::new(ErrorKind::Eof, "end of file"));
        }
        Ok(num_bytes)
    }
}

/// Retrieve name from the description line
fn get_name(description: &[u8]) -> Result<String> {
    if !is_description(description) {
        return Err(Error::new(ErrorKind::Input, "invalid fasta format"));
    }
    let description = String::from_utf8(description[1..].to_vec())?;
    Ok(description
        .trim_start_matches(SPACE)
        .split(SPACE)
        .take(1)
        .collect::<String>()
        .trim()
        .into())
}

/// Check to see if the line is a description line
fn is_description(line: &[u8]) -> bool {
    if line.is_empty() {
        return false;
    }
    if line[0] == SEQ_START {
        return true;
    }
    false
}

/// Count the number of bases in a line
fn count_bases(line: &[u8]) -> Result<usize> {
    Ok(std::str::from_utf8(line)?.trim().len())
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
            "{}",
            "Should work for example in documentation [Issue #14]",
        );
        assert_eq!(
            expected, record,
            "{}",
            "Should work for example in documentation [Issue #14]",
        );

        record.clear();
        let expected = Fai {
            name: "two".into(),
            length: 28,
            offset: 98,
            line_bases: 14,
            line_width: 15,
            qual_offset: None,
        };
        assert_eq!(
            reader.read(&mut record).unwrap_err().kind,
            ErrorKind::Eof,
            "{}",
            "Should return an end of file error [Issue #14]",
        );
        assert_eq!(
            expected, record,
            "{}",
            "Should work for example in documentation [Issue #14]",
        );
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
            "{}",
            "Should read a fasta entry into Fai record [Issue #14]"
        );
        assert_eq!(
            expected, record,
            "{}",
            "Should read in a record [Issue #14]",
        );

        record.clear();
        let expected = Fai {
            name: "abcdef".into(),
            offset: 25,
            line_width: 6,
            line_bases: 5,
            length: 5,
            qual_offset: None,
        };
        assert_eq!(
            ErrorKind::Eof,
            reader.read(&mut record).unwrap_err().kind,
            "{}",
            "Should return an Eof error [Issue #14]",
        );
        assert_eq!(
            expected, record,
            "{}",
            "Should read in a second record [Issue #14]",
        );
    }

    #[test]
    fn test_reader_read_description() {
        let input_line = std::io::Cursor::new(b">abcdef aaaa\nAAAA\nAAAA\n");
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        assert!(
            reader.read_description(&mut record).is_ok(),
            "{}",
            "Should read description into Fai record [Issue #14]",
        );
        assert_eq!(
            "abcdef", &record.name,
            "{}",
            "Should read in the record name [Issue #14]",
        );
        assert_eq!(
            13, record.offset,
            "{}",
            "Should read in the sequence offset [Issue #14]",
        );
    }

    #[test]
    fn test_reader_read_sequence() {
        let input_line = std::io::Cursor::new(b"AAAA\nAAAA\n>abcdef\nAAAA\n");
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        assert!(
            reader.read_sequence(&mut record).is_ok(),
            "{}",
            "Should read sequence into Fai record [Issue #14]"
        );
        assert_eq!(
            4, record.line_bases,
            "{}",
            "Should read in the number of bases in a line [Issue #14]",
        );
        assert_eq!(
            5, record.line_width,
            "{}",
            "Should read in the number of bytes in a line [Issue #14]",
        );
        assert_eq!(
            8, record.length,
            "{}",
            "Should increment number of bases [Issue #14]",
        );

        let input_line = std::io::Cursor::new(b"AAAA\nAAAAA\n>abcdef\nAAAA\n");
        let mut reader = Reader::new(input_line);
        record.clear();
        assert!(
            reader.read_sequence(&mut record).is_err(),
            "{}",
            "Should error out if the line widths are inconsistent [Issue #14]"
        );
    }

    #[test]
    fn test_reader_read_sequence_line() {
        let input_line = std::io::Cursor::new(b"AAAA\r\nAAAA\r\n");
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        assert!(
            reader.read_sequence_line(&mut record).is_ok(),
            "{}",
            "Should read sequence line into Fai record [Issue #14]"
        );
        assert_eq!(
            4, record.line_bases,
            "{}",
            "Should read in the number of bases in a line [Issue #14]",
        );
        assert_eq!(
            6, record.line_width,
            "{}",
            "Should read in the number of bytes in a line [Issue #14]",
        );
        assert_eq!(
            4, record.length,
            "{}",
            "Should increment number of bases [Issue #14]",
        );
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
                name: "Should return an error on an invalid description [Issue #14]",
                description: b"abcdefg",
                expect_error: true,
                expected: String::new(),
            },
            TestCase {
                name: "Should parse a name from the description [Issue #14]",
                description: b">name a b c d e f",
                expect_error: false,
                expected: "name".to_string(),
            },
            TestCase {
                name: "Should skip spaces after > before name [Issue #14]",
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
                name: "Should return true for a description line [Issue #14]",
                line: b">abcdef",
                expected: true,
            },
            TestCase {
                name: "Should return false for a non description line [Issue #14]",
                line: b"bla",
                expected: false,
            },
            TestCase {
                name: "Should return false for an empty line [Issue #14]",
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
    fn test_count_bases() {
        struct TestCase<'a> {
            name: &'a str,
            line: &'a [u8],
            expect_err: bool,
            expected: Result<usize>,
        }
        let test_cases = [TestCase {
            name: "Should return the base count [Issue #14]",
            line: b"AAAA\r\n",
            expect_err: false,
            expected: Ok(4),
        }];
        for test_case in test_cases {
            let actual = count_bases(test_case.line);
            if test_case.expect_err {
                assert!(actual.is_err(), "{}", test_case.name);
            } else {
                assert_eq!(test_case.expected, actual, "{}", test_case.name);
            }
        }
    }
}
