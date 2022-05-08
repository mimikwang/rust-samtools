use super::{Fai, IterFai, ReadFai, DELIMITER};
use crate::errors::{Error, ErrorKind, Result};

/// Reader is a reader for fai files
pub struct Reader<R: std::io::Read> {
    reader: csv::Reader<R>,
    string_record: csv::StringRecord,
}

impl<R> Reader<R>
where
    R: std::io::Read,
{
    /// Construct a fai reader from std::io::Read
    pub fn from_reader(reader: R) -> Self {
        Self {
            reader: csv::ReaderBuilder::new()
                .delimiter(DELIMITER)
                .has_headers(false)
                .from_reader(reader),
            string_record: csv::StringRecord::new(),
        }
    }

    // Returns an iterator
    pub fn iter(self) -> IterFai<Self> {
        IterFai::new(self)
    }
}

impl<R> ReadFai for Reader<R>
where
    R: std::io::Read,
{
    /// Read an fai record
    fn read(&mut self, record: &mut Fai) -> Result<()> {
        if self.reader.read_record(&mut self.string_record)? {
            *record = Fai::try_from(&mut self.string_record)?;
            return Ok(());
        }
        Err(Error::new(ErrorKind::Eof, "end of file"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reader_read() {
        struct TestCase<'a> {
            name: &'a str,
            input_line: &'a [u8],
            expect_error: bool,
            expected: Fai,
        }
        let test_cases = &mut [
            TestCase {
                name: "Should read in a FAI record [Issue #5]",
                input_line: b"name\t1\t2\t3\t4\n",
                expect_error: false,
                expected: Fai {
                    name: "name".into(),
                    length: 1,
                    offset: 2,
                    line_bases: 3,
                    line_width: 4,
                    qual_offset: None,
                },
            },
            TestCase {
                name: "Should error out if end of file [Issue #5]",
                input_line: b"",
                expect_error: true,
                expected: Fai::new(),
            },
            TestCase {
                name: "Should error out if the input is invalid [Issue #5]",
                input_line: b"name\tasdf\t2\t3\t4\t5\n",
                expect_error: true,
                expected: Fai::new(),
            },
        ];
        for test_case in test_cases.iter_mut() {
            let mut reader = Reader::from_reader(test_case.input_line);
            let mut record = Fai::new();
            let actual = reader.read(&mut record);
            if test_case.expect_error {
                assert!(actual.is_err(), "{}", test_case.name);
            } else {
                assert!(actual.is_ok(), "{}", test_case.name);
                assert_eq!(test_case.expected, record, "{}", test_case.name);
            }
        }
    }

    #[test]
    fn test_reader_read_different_widths() {
        let name = "Should error out if mixing fasta and fastq entries [Issue #5]";
        let input_lines: &[u8] = b"ref1\t1\t2\t3\t4\nref2\t1\t2\t3\t4\t5\n";
        let mut reader = Reader::from_reader(input_lines);
        let mut record = Fai::new();
        assert!(reader.read(&mut record).is_ok(), "{}", name);
        assert!(reader.read(&mut record).is_err(), "{}", name);

        let input_lines: &[u8] = b"ref1\t1\t2\t3\t4\t5\nref2\t1\t2\t3\t4\n";
        let mut reader = Reader::from_reader(input_lines);
        let mut record = Fai::new();
        assert!(reader.read(&mut record).is_ok(), "{}", name);
        assert!(reader.read(&mut record).is_err(), "{}", name);
    }

    #[test]
    fn test_reader_iter() {
        let name = "Should iterate through multiple lines [Issue #5]";
        let input_lines: &[u8] = b"ref1\t1\t2\t3\t4\t5\nref2\t1\t2\t3\t4\t5\n";
        let reader = Reader::from_reader(input_lines);
        let mut iter = reader.iter();

        assert_eq!(
            iter.next(),
            Some(Ok(Fai {
                name: "ref1".into(),
                length: 1,
                offset: 2,
                line_bases: 3,
                line_width: 4,
                qual_offset: Some(5),
            })),
            "{}",
            name,
        );
        assert_eq!(
            iter.next(),
            Some(Ok(Fai {
                name: "ref2".into(),
                length: 1,
                offset: 2,
                line_bases: 3,
                line_width: 4,
                qual_offset: Some(5),
            })),
            "{}",
            name,
        );
        assert_eq!(iter.next(), None, "{}", name);
    }
}
