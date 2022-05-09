/// Shorthand for Result<T, Error>
pub type Result<T> = std::result::Result<T, Error>;

/// Custom error type for the crate
///
/// All external error types should be converted to `Error`.
///
#[derive(Debug, PartialEq)]
pub struct Error {
    pub kind: ErrorKind,
    pub message: String,
}

/// Kind of error
#[derive(Debug, PartialEq)]
pub enum ErrorKind {
    /// Input format errors
    Input,
    /// Errors related to input / output
    IO,
    /// End of file errors
    Eof,
    /// Type conversion errors
    TypeConversion,
    /// User related errors
    User,
    /// All other errors
    Unknown,
}

impl Default for ErrorKind {
    fn default() -> Self {
        Self::Unknown
    }
}

impl std::error::Error for Error {}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "kind: {:?}, message: {}", self.kind, self.message)
    }
}

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
