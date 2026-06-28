//! Memory-mapped flat fingerprint store for O(1) lookup by doc_id.
//!
//! The store is a contiguous array of fingerprints written in doc_id order —
//! the same order the Lance dataset is scanned when building the inverted index,
//! so `doc_id` indexes directly into the file. Each fingerprint occupies
//! `FP_WORDS * 8 = 128` bytes (16 little-endian u64 words).
//!
//! For the full Enamine REAL 1.4B set this is ~174 GB on disk; mmap keeps the
//! resident set bounded by the OS page cache while BMW touches only the small
//! fraction of fingerprints that survive pruning.

use std::path::Path;

use memmap2::Mmap;

use crate::error::{BitMakoError, Result};
use crate::etl::fingerprint::{Fingerprint, FP_WORDS};

/// Bytes per stored fingerprint (16 × u64).
pub const FP_BYTES: usize = FP_WORDS * 8;

/// Read-only memory-mapped view over a flat fingerprint file.
pub struct FpStore {
    mmap: Mmap,
    num_fps: usize,
}

impl FpStore {
    /// Open a flat fingerprint store produced by `build-fp-store`.
    pub fn open(path: &Path) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let len = mmap.len();
        if len % FP_BYTES != 0 {
            return Err(BitMakoError::IndexBuild(format!(
                "fingerprint store size {} is not a multiple of {} bytes",
                len, FP_BYTES
            )));
        }
        Ok(FpStore { mmap, num_fps: len / FP_BYTES })
    }

    /// Number of fingerprints in the store.
    #[inline]
    pub fn len(&self) -> usize {
        self.num_fps
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.num_fps == 0
    }

    /// Fetch the fingerprint for `doc_id`, or `None` if out of range.
    ///
    /// Bytes are copied out via `from_le_bytes`, so there are no alignment
    /// requirements on the underlying mmap.
    #[inline]
    pub fn get(&self, doc_id: u32) -> Option<Fingerprint> {
        let idx = doc_id as usize;
        if idx >= self.num_fps {
            return None;
        }
        let start = idx * FP_BYTES;
        let bytes = &self.mmap[start..start + FP_BYTES];
        let mut fp = [0u64; FP_WORDS];
        for (i, chunk) in bytes.chunks_exact(8).enumerate() {
            fp[i] = u64::from_le_bytes(chunk.try_into().unwrap());
        }
        Some(fp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::fingerprint::compute_morgan_fp;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_store(fps: &[Fingerprint]) -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        for fp in fps {
            let mut buf = [0u8; FP_BYTES];
            for (i, &w) in fp.iter().enumerate() {
                buf[i * 8..i * 8 + 8].copy_from_slice(&w.to_le_bytes());
            }
            tmp.write_all(&buf).unwrap();
        }
        tmp.flush().unwrap();
        tmp
    }

    #[test]
    fn test_roundtrip() {
        let fps = vec![
            compute_morgan_fp("CCO"),
            compute_morgan_fp("c1ccccc1"),
            compute_morgan_fp("CNC(=O)c1ccccc1"),
        ];
        let tmp = write_store(&fps);
        let store = FpStore::open(tmp.path()).unwrap();

        assert_eq!(store.len(), 3);
        for (i, fp) in fps.iter().enumerate() {
            assert_eq!(store.get(i as u32).as_ref(), Some(fp));
        }
        assert_eq!(store.get(3), None);
    }
}
