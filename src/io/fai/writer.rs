use super::{Fai, DELIMITER};
use crate::errors::Result;

/// Writer is a writer for fai files
pub struct Writer<W: std::io::Write> {
    writer: csv::Writer<W>,
}

impl<W> Writer<W>
where
    W: std::io::Write,
{
    /// Construct a Fai writer from `std::io::Write`
    pub fn new(writer: W) -> Self {
        Self {
            writer: csv::WriterBuilder::new()
                .delimiter(DELIMITER)
                .has_headers(false)
                .from_writer(writer),
        }
    }

    /// Write a Fai record
    pub fn write(&mut self, record: &Fai) -> Result<()> {
        self.writer.write_record(&record.to_string_record())?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_writer_write() {
        struct TestCase<'a> {
            name: &'a str,
            record: Fai,
            expected: &'a str,
        }
        let test_cases = &[
            TestCase {
                name: "Should write a fasta FAI record",
                record: Fai {
                    name: "record1".into(),
                    length: 1,
                    offset: 2,
                    line_bases: 3,
                    line_width: 4,
                    qual_offset: None,
                },
                expected: "record1\t1\t2\t3\t4\n",
            },
            TestCase {
                name: "Should write a fastq FAI record",
                record: Fai {
                    name: "record2".into(),
                    length: 1,
                    offset: 2,
                    line_bases: 3,
                    line_width: 4,
                    qual_offset: Some(5),
                },
                expected: "record2\t1\t2\t3\t4\t5\n",
            },
        ];
        for test_case in test_cases.iter() {
            let mut writer = Writer::new(vec![]);
            assert!(
                writer.write(&test_case.record).is_ok(),
                "{}",
                test_case.name
            );
            let data = String::from_utf8(writer.writer.into_inner().unwrap()).unwrap();
            assert_eq!(&data, test_case.expected, "{}", test_case.name);
        }
    }

    #[test]
    fn test_writer_write_different_widths() {
        let name = "Should error out if mixing fasta and fastq entries";
        let mut writer = Writer::new(vec![]);
        let mut record = Fai {
            name: "name".into(),
            length: 1,
            offset: 2,
            line_bases: 3,
            line_width: 4,
            qual_offset: None,
        };
        assert!(writer.write(&record).is_ok(), "{}", name);
        record.qual_offset = Some(1);
        assert!(writer.write(&record).is_err(), "{}", name);
    }
}
