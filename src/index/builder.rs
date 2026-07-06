//! Inverted index builder for 1024-bit fingerprint bit vectors.
//!
//! Processes parsed compounds in streaming fashion to build 1024 posting lists,
//! one per fingerprint bit. After all compounds are ingested, the index is
//! finalized (sorted, block-max populated) and serialized to disk.

use std::io::Write;
use std::path::Path;
use std::fs::File;

use tracing::info;

use crate::error::{BitMakoError, Result};
use crate::etl::fingerprint::{fp_popcount, Fingerprint, FP_BITS};
use crate::index::posting_list::{PostingList, BLOCK_SIZE};
use crate::index::{INDEX_MAGIC, INDEX_VERSION};

/// Accumulates compound fingerprints and builds posting lists incrementally.
pub struct IndexBuilder {
    /// posting_lists[bit] = sorted list of doc_ids with that bit set
    posting_lists: Vec<Vec<u32>>,
    /// Popcount for each compound (indexed by doc_id)
    compound_pops: Vec<u8>,
    /// Total compounds ingested
    num_compounds: u32,
}

impl Default for IndexBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl IndexBuilder {
    pub fn new() -> Self {
        IndexBuilder {
            posting_lists: vec![Vec::new(); FP_BITS],
            compound_pops: Vec::new(),
            num_compounds: 0,
        }
    }

    /// Add a compound fingerprint to the index.
    /// `doc_id` must be monotonically increasing (for sorted posting lists).
    pub fn add_compound(&mut self, doc_id: u32, fingerprint: &Fingerprint) {
        let pop = fp_popcount(fingerprint) as u8;
        self.compound_pops.push(pop);

        // Set bits in the 1024-bit fingerprint, push doc_id to each posting list
        for (word_idx, &word) in fingerprint.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let trailing = bits.trailing_zeros() as usize;
                let bit_idx = word_idx * 64 + trailing;
                self.posting_lists[bit_idx].push(doc_id);
                bits &= bits - 1; // clear lowest set bit
            }
        }

        self.num_compounds += 1;

        if self.num_compounds % 1_000_000 == 0 {
            info!("Index builder: {} million compounds ingested", self.num_compounds / 1_000_000);
        }
    }

    /// Finalize posting lists (they should already be sorted since doc_ids
    /// are monotonically increasing, but we sort just in case).
    fn finalize_lists(&mut self) -> Vec<PostingList> {
        let compound_pops = &self.compound_pops;
        let mut finalized = Vec::with_capacity(FP_BITS);

        for ids in self.posting_lists.iter_mut() {
            ids.sort_unstable();

            let num_blocks = ids.len().div_ceil(BLOCK_SIZE);
            let mut block_max_pop: Vec<u8> = Vec::with_capacity(num_blocks);

            for block_idx in 0..num_blocks {
                let start = block_idx * BLOCK_SIZE;
                let end = (start + BLOCK_SIZE).min(ids.len());
                let max_pop = ids[start..end]
                    .iter()
                    .map(|&doc_id| compound_pops.get(doc_id as usize).copied().unwrap_or(0))
                    .max()
                    .unwrap_or(0);
                block_max_pop.push(max_pop);
            }

            finalized.push(PostingList {
                doc_ids: std::mem::take(ids),
                block_max_pop,
            });
        }

        finalized
    }

    /// Build and write the finalized index to disk.
    ///
    /// File layout:
    ///   [8-byte magic]
    ///   [4-byte version]
    ///   [4-byte num_compounds]
    ///   [4-byte num_bits = 1024]
    ///   [8-byte * 1024: byte offsets to each posting list from start of postings section]
    ///   [4-byte * num_compounds: compound_pops, packed as u8 → padded to u32 alignment]
    ///   [posting list data for bits 0..1023]
    pub fn write_index(mut self, output_path: &Path) -> Result<IndexStats> {
        info!("Finalizing index for {} compounds", self.num_compounds);
        let lists = self.finalize_lists();

        let mut file = File::create(output_path)?;

        // Serialize all posting lists
        let mut serialized_lists: Vec<Vec<u8>> = Vec::with_capacity(FP_BITS);
        for pl in &lists {
            serialized_lists.push(pl.serialize().map_err(BitMakoError::Io)?);
        }

        // Compute offsets (u64 to support >4 GiB posting sections)
        let mut offsets: Vec<u64> = Vec::with_capacity(FP_BITS);
        let mut running_offset: u64 = 0;
        for data in &serialized_lists {
            offsets.push(running_offset);
            running_offset += data.len() as u64;
        }

        // Write header
        file.write_all(INDEX_MAGIC)?;
        file.write_all(&INDEX_VERSION.to_le_bytes())?;
        file.write_all(&self.num_compounds.to_le_bytes())?;
        file.write_all(&(FP_BITS as u32).to_le_bytes())?;

        // Write offsets table (8 bytes each)
        for off in &offsets {
            file.write_all(&off.to_le_bytes())?;
        }

        // Write compound popcounts (one u8 per compound)
        file.write_all(&self.compound_pops)?;
        // Pad to 4-byte alignment
        let pad = (4 - (self.compound_pops.len() % 4)) % 4;
        for _ in 0..pad {
            file.write_all(&[0u8])?;
        }

        // Write posting lists
        let mut total_postings = 0usize;
        for data in &serialized_lists {
            file.write_all(data)?;
            total_postings += data.len();
        }

        let stats = IndexStats {
            num_compounds: self.num_compounds as usize,
            total_postings_bytes: total_postings,
            non_empty_bits: lists.iter().filter(|pl| !pl.is_empty()).count(),
        };

        info!(
            "Index written: compounds={} bits_active={} postings_bytes={}",
            stats.num_compounds, stats.non_empty_bits, stats.total_postings_bytes
        );

        Ok(stats)
    }
}

/// Summary statistics from an index build
#[derive(Debug)]
pub struct IndexStats {
    pub num_compounds: usize,
    pub total_postings_bytes: usize,
    pub non_empty_bits: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::fingerprint::compute_morgan_fp;
    use tempfile::NamedTempFile;

    #[test]
    fn test_builder_add_and_write() {
        let mut builder = IndexBuilder::new();
        let fp1 = compute_morgan_fp("CCO");
        let fp2 = compute_morgan_fp("c1ccccc1");
        builder.add_compound(0, &fp1);
        builder.add_compound(1, &fp2);

        let tmp = NamedTempFile::new().unwrap();
        let stats = builder.write_index(tmp.path()).unwrap();
        assert_eq!(stats.num_compounds, 2);
        assert!(stats.non_empty_bits > 0);
        assert!(stats.total_postings_bytes > 0);
    }
}
