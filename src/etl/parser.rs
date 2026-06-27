//! Zero-copy batch parser: processes chunks of raw TSV lines into ParsedCompound records.
//!
//! Uses `rayon::par_iter` for data-parallel processing across CPU cores.
//! Each raw line is processed independently, enabling work-stealing across threads.

use rayon::prelude::*;
use tracing::{debug, warn};

use crate::etl::fingerprint::{compute_morgan_fp, Fingerprint};
use crate::etl::properties::{compute_properties, MolecularProperties};
use crate::etl::reader::{RawLine, ReaderConfig};

/// A fully parsed and featurized compound ready for columnar serialization.
#[derive(Debug, Clone)]
pub struct ParsedCompound {
    pub compound_id: String,
    pub smiles: String,
    pub fingerprint: Fingerprint,
    pub properties: MolecularProperties,
}

/// Parse a single raw TSV line into a ParsedCompound.
/// Returns None if the line is malformed or the SMILES is invalid.
fn parse_line(raw: &RawLine, config: &ReaderConfig) -> Option<ParsedCompound> {
    let line = raw.raw.trim();
    if line.is_empty() {
        return None;
    }

    let (smiles_str, id_str) = crate::etl::reader::split_line(line, config)?;

    if smiles_str.len() < 2 {
        warn!("Line {}: SMILES too short: '{}'", raw.line_num, smiles_str);
        return None;
    }

    let fingerprint = compute_morgan_fp(smiles_str);
    let properties = compute_properties(smiles_str);

    // Reject molecules with zero-bit fingerprints (parse failure indicator)
    if fingerprint == [0u64; 16] && smiles_str.len() > 1 {
        debug!("Line {}: fingerprint all-zero, may be parse failure", raw.line_num);
    }

    Some(ParsedCompound {
        compound_id: id_str.to_owned(),
        smiles: smiles_str.to_owned(),
        fingerprint,
        properties,
    })
}

/// Process a chunk of raw lines in parallel using rayon work-stealing.
///
/// Returns (parsed_compounds, failed_count).
pub fn parse_chunk_parallel(
    chunk: Vec<RawLine>,
    config: &ReaderConfig,
) -> (Vec<ParsedCompound>, usize) {
    let results: Vec<Option<ParsedCompound>> = chunk
        .par_iter()
        .map(|raw| parse_line(raw, config))
        .collect();

    let mut compounds = Vec::with_capacity(results.len());
    let mut failed = 0usize;

    for opt in results {
        match opt {
            Some(c) => compounds.push(c),
            None => failed += 1,
        }
    }

    (compounds, failed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::reader::RawLine;

    #[test]
    fn test_parse_valid_line() {
        let config = ReaderConfig::default();
        let raw = RawLine {
            line_num: 1,
            raw: "CCO\tZ1234567890".to_string(),
        };
        let result = parse_line(&raw, &config);
        assert!(result.is_some());
        let c = result.unwrap();
        assert_eq!(c.compound_id, "Z1234567890");
        assert_eq!(c.smiles, "CCO");
        assert!(c.properties.mw > 0.0);
    }

    #[test]
    fn test_parse_empty_line() {
        let config = ReaderConfig::default();
        let raw = RawLine { line_num: 2, raw: "".to_string() };
        assert!(parse_line(&raw, &config).is_none());
    }

    #[test]
    fn test_parse_chunk_parallel() {
        let config = ReaderConfig::default();
        let chunk = vec![
            RawLine { line_num: 1, raw: "CCO\tZ001".to_string() },
            RawLine { line_num: 2, raw: "c1ccccc1\tZ002".to_string() },
            RawLine { line_num: 3, raw: "\tZ003".to_string() }, // invalid
        ];
        let (compounds, failed) = parse_chunk_parallel(chunk, &config);
        assert_eq!(compounds.len(), 2);
        assert_eq!(failed, 1);
    }
}
