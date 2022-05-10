use super::{ReadToFai, Record, Records, DELIMITER};
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
    /// Construct a Fai reader from `std::io::Read`
    pub fn new(reader: R) -> Self {
        Self {
            reader: csv::ReaderBuilder::new()
                .delimiter(DELIMITER)
                .has_headers(false)
                .from_reader(reader),
            string_record: csv::StringRecord::new(),
        }
    }

    // Consume the reader and return an `IterFai`
    pub fn iter(self) -> Records<Self> {
        Records::new(self)
    }
}

impl<R> ReadToFai for Reader<R>
where
    R: std::io::Read,
{
    /// Read a Fai record
    fn read(&mut self, record: &mut Record) -> Result<()> {
        if self.reader.read_record(&mut self.string_record)? {
            *record = Record::try_from(&mut self.string_record)?;
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
            expected: Record,
        }
        let test_cases = &mut [
            TestCase {
                name: "Should read in a FAI record",
                input_line: b"name\t1\t2\t3\t4\n",
                expect_error: false,
                expected: Record {
                    name: "name".into(),
                    length: 1,
                    offset: 2,
                    line_bases: 3,
                    line_width: 4,
                    qual_offset: None,
                },
            },
            TestCase {
                name: "Should error out if end of file",
                input_line: b"",
                expect_error: true,
                expected: Record::new(),
            },
            TestCase {
                name: "Should error out if the input is invalid",
                input_line: b"name\tasdf\t2\t3\t4\t5\n",
                expect_error: true,
                expected: Record::new(),
            },
        ];
        for test_case in test_cases.iter_mut() {
            let mut reader = Reader::new(test_case.input_line);
            let mut record = Record::new();
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
        let name = "Should error out if mixing fasta and fastq entries";
        let input_lines: &[u8] = b"ref1\t1\t2\t3\t4\nref2\t1\t2\t3\t4\t5\n";
        let mut reader = Reader::new(input_lines);
        let mut record = Record::new();
        assert!(reader.read(&mut record).is_ok(), "{}", name);
        assert!(reader.read(&mut record).is_err(), "{}", name);

        let input_lines: &[u8] = b"ref1\t1\t2\t3\t4\t5\nref2\t1\t2\t3\t4\n";
        let mut reader = Reader::new(input_lines);
        let mut record = Record::new();
        assert!(reader.read(&mut record).is_ok(), "{}", name);
        assert!(reader.read(&mut record).is_err(), "{}", name);
    }

    #[test]
    fn test_reader_iter() {
        let input_lines: &[u8] = b"ref1\t1\t2\t3\t4\t5\nref2\t1\t2\t3\t4\t5\n";
        let reader = Reader::new(input_lines);
        let mut iter = reader.iter();

        assert_eq!(
            iter.next(),
            Some(Ok(Record {
                name: "ref1".into(),
                length: 1,
                offset: 2,
                line_bases: 3,
                line_width: 4,
                qual_offset: Some(5),
            })),
            "Should iterate and return the first entry",
        );
        assert_eq!(
            iter.next(),
            Some(Ok(Record {
                name: "ref2".into(),
                length: 1,
                offset: 2,
                line_bases: 3,
                line_width: 4,
                qual_offset: Some(5),
            })),
            "Should iterate and return the second entry",
        );
        assert_eq!(
            iter.next(),
            None,
            "Should return None if the iterator is consumed",
        );
    }
}
