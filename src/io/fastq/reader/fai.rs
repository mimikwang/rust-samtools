use super::super::super::common;
use super::*;
use crate::errors::{Error, ErrorKind, Result};
use crate::io::fai::{Fai, IterFai, ReadFai};

/// Reader reads a fastq file into FAI entries
pub struct Reader<B>
where
    B: std::io::BufRead + std::io::Seek,
{
    reader: B,
    buffer: Vec<u8>,
    sequence_num_bytes: usize,
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
            sequence_num_bytes: 0,
        }
    }

    // Returns an iterator
    pub fn iter(self) -> IterFai<Self> {
        IterFai::new(self)
    }

    /// Read the first line of the fastq entry
    fn read_description(&mut self, record: &mut Fai) -> Result<()> {
        if self.buffer.is_empty() {
            self.read_line()?;
        }
        record.name = get_name(&self.buffer)?;
        record.offset = self.reader.stream_position()?;
        self.buffer.clear();
        Ok(())
    }

    /// Read the entire sequence
    fn read_sequence(&mut self, record: &mut Fai) -> Result<()> {
        self.sequence_num_bytes = 0;
        loop {
            if is_plus(&self.buffer) {
                return Ok(());
            }
            self.read_sequence_line(record)?;
        }
    }

    /// Read in a sequence line
    fn read_sequence_line(&mut self, record: &mut Fai) -> Result<()> {
        self.buffer.clear();
        let num_bytes = self.read_line()?;
        if is_plus(&self.buffer) {
            return Ok(());
        }
        if record.line_width == 0 {
            record.line_width = num_bytes;
            record.line_bases = common::count_bases(&self.buffer)?;
        } else if record.line_width < num_bytes {
            return Err(Error::new(ErrorKind::Input, "invalid fastq record"));
        }
        self.sequence_num_bytes += num_bytes;
        record.length += common::count_bases(&self.buffer)?;
        Ok(())
    }

    /// Read in the plus line
    fn read_plus(&mut self, record: &mut Fai) -> Result<()> {
        record.qual_offset = Some(self.reader.stream_position()?);
        self.buffer.clear();
        Ok(())
    }

    /// Read the quality portion
    fn read_quality(&mut self) -> Result<()> {
        self.reader.consume(self.sequence_num_bytes);
        self.read_line()?;
        Ok(())
    }

    /// Read in a line of data
    fn read_line(&mut self) -> Result<usize> {
        common::read_line(&mut self.reader, &mut self.buffer)
    }
}

impl<R> ReadFai for Reader<R>
where
    R: std::io::BufRead + std::io::Seek,
{
    /// Read a FAI record
    fn read(&mut self, record: &mut Fai) -> Result<()> {
        self.read_description(record)?;
        self.read_sequence(record)?;
        self.read_plus(record)?;
        self.read_quality()?;
        Ok(())
    }
}

/// Retrieve name from the description line
fn get_name(description: &[u8]) -> Result<String> {
    if !is_description(description) {
        return Err(Error::new(ErrorKind::Input, "invalid fastq format"));
    }
    let description = String::from_utf8(description[1..].to_vec())?;
    Ok(common::parse_sequence_name(&description))
}

/// Check to see if the line is a description line
fn is_description(line: &[u8]) -> bool {
    line.starts_with(&[SEQ_START])
}

/// Check to see if the line is a plus line
fn is_plus(line: &[u8]) -> bool {
    line.starts_with(&[PLUS])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_example_from_samtools_doc() {
        // Test example from samtools documentation here:
        //  https://www.htslib.org/doc/faidx.html
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
        let mut reader = Reader::new(input_line);
        let mut record = Fai::new();
        let expected = Fai {
            name: "fastq1".into(),
            length: 66,
            offset: 8,
            line_bases: 30,
            line_width: 31,
            qual_offset: Some(79),
        };
        assert!(
            reader.read(&mut record).is_ok(),
            "{}",
            "Should work for example in documentation [Issue #15]",
        );
        assert_eq!(
            expected, record,
            "{}",
            "Should work for example in documentation [Issue #15]",
        );

        record.clear();
        let expected = Fai {
            name: "fastq2".into(),
            length: 28,
            offset: 156,
            line_bases: 14,
            line_width: 15,
            qual_offset: Some(188),
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
    fn test_is_description() {
        struct TestCase<'a> {
            name: &'a str,
            line: &'a [u8],
            expected: bool,
        }
        let test_cases = [
            TestCase {
                name: "Should return true for a description line [Issue #15]",
                line: b"@abcdef",
                expected: true,
            },
            TestCase {
                name: "Should return false for a non description line [Issue #15]",
                line: b"bla",
                expected: false,
            },
            TestCase {
                name: "Should return false for an empty line [Issue #15]",
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
}
