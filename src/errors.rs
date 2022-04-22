/// Shorthand for Result<T, Error>
type Result<T> = std::result::Result<T, Error>;

/// Custom error type for the crate
#[derive(Debug, PartialEq)]
pub struct Error {
    pub kind: ErrorKind,
    pub message: String,
}

/// Kind of error
#[derive(Debug, PartialEq)]
pub enum ErrorKind {
    UnknownError,
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
