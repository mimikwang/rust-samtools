/// Shorthand for Result<T, Error>
pub type Result<T> = std::result::Result<T, Error>;

/// Custom error type for the crate
#[derive(Debug, PartialEq)]
pub struct Error {
    pub kind: ErrorKind,
    pub message: String,
}

/// Kind of error
#[derive(Debug, PartialEq)]
pub enum ErrorKind {
    Input,
    IO,
    Eof,
    TypeConversion,
    Unknown,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "kind: {:?}, message: {}", self.kind, self.message)
    }
}

impl std::error::Error for Error {}

impl Error {
    /// Basic constructor for Error
    pub fn new(kind: ErrorKind, message: &str) -> Self {
        Error {
            kind,
            message: message.into(),
        }
    }
}

impl From<csv::Error> for Error {
    fn from(e: csv::Error) -> Self {
        Self::new(ErrorKind::Input, &e.to_string())
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Self::new(ErrorKind::IO, &e.to_string())
    }
}

impl From<std::string::FromUtf8Error> for Error {
    fn from(e: std::string::FromUtf8Error) -> Self {
        Self::new(ErrorKind::TypeConversion, &e.to_string())
    }
}

impl From<std::str::Utf8Error> for Error {
    fn from(e: std::str::Utf8Error) -> Self {
        Self::new(ErrorKind::TypeConversion, &e.to_string())
    }
}
