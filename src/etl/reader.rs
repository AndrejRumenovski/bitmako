//! Streaming bzip2 reader that decompresses on-the-fly and yields compound chunks.
//!
//! Never reads the full file into RAM. Chunks of lines are accumulated then
//! sent to the parallel processing pipeline via a crossbeam channel.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use bzip2::read::BzDecoder;
use crossbeam_channel::Sender;
use tracing::{debug, info};

use crate::error::{BitMakoError, Result};

/// Lines-of-text owned batch ready for parallel processing
pub type LineChunk = Vec<String>;

/// Configuration for the streaming reader
#[derive(Debug, Clone)]
pub struct ReaderConfig {
    /// Number of compound lines per processing chunk
    pub chunk_size: usize,
    /// Maximum number of in-flight chunks in the channel backpressure buffer
    pub channel_capacity: usize,
    /// Tab-separated column index for SMILES (default 0 for Enamine REAL format)
    pub smiles_col: usize,
    /// Tab-separated column index for compound ID (default 1)
    pub id_col: usize,
    /// Skip header line
    pub has_header: bool,
}

impl Default for ReaderConfig {
    fn default() -> Self {
        ReaderConfig {
            chunk_size: 100_000,
            channel_capacity: 8,
            smiles_col: 0,
            id_col: 1,
            has_header: true,
        }
    }
}

/// A raw compound line before parsing
#[derive(Debug, Clone)]
pub struct RawLine {
    pub line_num: usize,
    pub raw: String,
}

/// Stream lines from a `.bz2` compressed TSV file.
/// Opens the file, wraps it in a BzDecoder, and sends chunks of raw lines
/// into the provided crossbeam channel for parallel downstream processing.
///
/// Returns the total number of lines sent (excluding header).
pub fn stream_bz2_file(
    path: &Path,
    config: &ReaderConfig,
    tx: Sender<Vec<RawLine>>,
) -> Result<usize> {
    let file = File::open(path)?;
    let decoder = BzDecoder::new(file);
    let reader = BufReader::with_capacity(4 * 1024 * 1024, decoder); // 4 MB read buffer

    let mut lines = reader.lines();
    let mut line_num = 0usize;
    let mut total_sent = 0usize;
    let mut chunk: Vec<RawLine> = Vec::with_capacity(config.chunk_size);

    // Skip header if present
    if config.has_header {
        if let Some(header) = lines.next() {
            let _ = header?;
            line_num += 1;
            debug!("Skipped header line");
        }
    }

    for line_result in lines {
        let line = line_result?;
        line_num += 1;

        if line.trim().is_empty() {
            continue;
        }

        chunk.push(RawLine { line_num, raw: line });

        if chunk.len() >= config.chunk_size {
            let batch = std::mem::replace(&mut chunk, Vec::with_capacity(config.chunk_size));
            let batch_size = batch.len();
            tx.send(batch).map_err(|e| BitMakoError::ChannelSend(e.to_string()))?;
            total_sent += batch_size;
            info!("Streamed {} compounds (total {})", batch_size, total_sent);
        }
    }

    // Send final partial chunk
    if !chunk.is_empty() {
        total_sent += chunk.len();
        tx.send(chunk).map_err(|e| BitMakoError::ChannelSend(e.to_string()))?;
    }

    info!("Streaming complete. Total lines: {}", total_sent);
    Ok(total_sent)
}

/// Extract SMILES and compound ID from a raw TSV line using column config.
/// Returns borrowed slices into the line to avoid allocation.
pub fn split_line<'a>(line: &'a str, config: &ReaderConfig) -> Option<(&'a str, &'a str)> {
    let mut cols = line.splitn(3, '\t');
    let smiles = match config.smiles_col {
        0 => cols.next()?.trim(),
        1 => { cols.next(); cols.next()?.trim() }
        _ => return None,
    };

    // Re-split cleanly for the ID column
    let mut cols2 = line.splitn(3, '\t');
    let id = match config.id_col {
        0 => cols2.next()?.trim(),
        1 => { cols2.next(); cols2.next()?.trim() }
        2 => { cols2.next(); cols2.next(); cols2.next()?.trim() }
        _ => return None,
    };

    if smiles.is_empty() {
        return None;
    }

    Some((smiles, id))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_line() {
        let cfg = ReaderConfig::default();
        let line = "CCO\tZ1234567890\t";
        let (smiles, id) = split_line(line, &cfg).unwrap();
        assert_eq!(smiles, "CCO");
        assert_eq!(id, "Z1234567890");
    }

    #[test]
    fn test_split_empty_smiles() {
        let cfg = ReaderConfig::default();
        let line = "\tZ1234567890";
        assert!(split_line(line, &cfg).is_none());
    }
}
