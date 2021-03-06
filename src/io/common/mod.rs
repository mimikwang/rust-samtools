use crate::errors::{Error, ErrorKind, Result};

const SPACE: char = ' ';
const NEWLINE: u8 = b'\n';

/// Parse sequence name from a line
///
/// A sequence name is defined as the part of a line up to the first space.
///
pub fn parse_sequence_name(line: &str) -> String {
    line.trim()
        .split(SPACE)
        .take(1)
        .collect::<String>()
        .trim()
        .into()
}

/// Count the number of bases in a line
///
/// The number of bases in a line is defined as any byte that isn't a leading or trailing white
/// space.
///
pub fn count_bases(line: &[u8]) -> Result<usize> {
    Ok(std::str::from_utf8(line)?.trim().len())
}

/// Read a line of data
///
/// Read a line of data into `buffer` and return the number of bytes read.  An end of file error
/// is returned if the number of bytes read is 0.
///
pub fn read_line<B>(reader: &mut B, buffer: &mut Vec<u8>) -> Result<usize>
where
    B: std::io::BufRead,
{
    let num_bytes = reader.read_until(NEWLINE, buffer)?;
    if num_bytes == 0 {
        return Err(Error::new(ErrorKind::Eof, "end of file"));
    }
    Ok(num_bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sequence_name() {
        struct TestCase<'a> {
            name: &'a str,
            line: &'a str,
            expected: &'a str,
        }
        let test_cases = [
            TestCase {
                name: "Should parse a name from a line",
                line: "name description a b c",
                expected: "name",
            },
            TestCase {
                name: "Should remove trailing new line",
                line: "name\n",
                expected: "name",
            },
            TestCase {
                name: "Should remove leading spaces",
                line: "  abc",
                expected: "abc",
            },
        ];
        for test_case in test_cases {
            assert_eq!(
                String::from(test_case.expected),
                parse_sequence_name(test_case.line),
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
            name: "Should return the base count",
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

    #[test]
    fn test_read_line() {
        let mut input: &[u8] = b"abcdefg";
        let mut buffer = Vec::new();
        assert_eq!(
            Ok(7),
            read_line(&mut input, &mut buffer),
            "Should read the correct number of bytes",
        );
        assert_eq!(
            ErrorKind::Eof,
            read_line(&mut input, &mut buffer).unwrap_err().kind,
            "Should return an Eof error",
        );
    }
}
