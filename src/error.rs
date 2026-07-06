use thiserror::Error;

#[derive(Error, Debug)]
pub enum BitMakoError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("Invalid SMILES at line {line}: '{smiles}'")]
    InvalidSmiles { line: usize, smiles: String },

    #[error("Arrow error: {0}")]
    Arrow(#[from] arrow_schema::ArrowError),

    #[error("Lance error: {0}")]
    Lance(String),

    #[error("Index build error: {0}")]
    IndexBuild(String),

    #[error("Index corrupt: expected block count {expected}, got {actual}")]
    IndexCorrupt { expected: usize, actual: usize },

    #[error("Query error: {0}")]
    Query(String),

    #[error("Serialization error: {0}")]
    Serialization(String),

    #[error("Channel send error: {0}")]
    ChannelSend(String),
}

impl From<bincode::Error> for BitMakoError {
    fn from(e: bincode::Error) -> Self {
        BitMakoError::Serialization(e.to_string())
    }
}

pub type Result<T> = std::result::Result<T, BitMakoError>;

/// Converts any displayable error (Lance's error type in practice) into
/// `BitMakoError::Lance`. Lance operations (`Dataset::open`, `.scan()`, `.take()`,
/// stream iteration, ...) all return their own error type, so without this every
/// call site needs `.map_err(|e| BitMakoError::Lance(e.to_string()))` spelled out —
/// this lets call sites write `.lance_err()?` instead.
pub trait LanceResultExt<T> {
    fn lance_err(self) -> Result<T>;
}

impl<T, E: std::fmt::Display> LanceResultExt<T> for std::result::Result<T, E> {
    fn lance_err(self) -> Result<T> {
        self.map_err(|e| BitMakoError::Lance(e.to_string()))
    }
}
