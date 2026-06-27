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
