//! Memory-mapped flat molecular-property store for O(1) lookup by doc_id.
//!
//! Layout mirrors the flat fingerprint store: properties are written in doc_id
//! order — the same order the Lance dataset is scanned when building the inverted
//! index and the fingerprint store — so `doc_id` indexes directly into the file.
//! Each record holds five little-endian `f32` values (`mw`, `logp`, `rot_bonds`,
//! `heavy_atoms`, `ring_count`) for a fixed `PROP_BYTES = 20` bytes per compound.
//!
//! The integer-valued fields are stored as `f32` because they are small enough to
//! round-trip exactly and `PropertyFilter` compares them as `f32` anyway. For the
//! full Enamine REAL 1.4B set this is ~27 GB on disk — small enough that the OS
//! page cache keeps the working set hot, letting BMW screen properties cheaply
//! before paying for the 128-byte fingerprint fetch.

use std::path::Path;

use memmap2::Mmap;

use crate::error::{BitMakoError, Result};
use crate::etl::properties::MolecularProperties;

/// Bytes per stored property record (5 × f32).
pub const PROP_BYTES: usize = 5 * 4;

/// Read-only memory-mapped view over a flat property file.
pub struct PropStore {
    mmap: Mmap,
    num_records: usize,
}

impl PropStore {
    /// Open a flat property store produced by `build-prop-store`.
    pub fn open(path: &Path) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let len = mmap.len();
        if len % PROP_BYTES != 0 {
            return Err(BitMakoError::IndexBuild(format!(
                "property store size {} is not a multiple of {} bytes",
                len, PROP_BYTES
            )));
        }
        Ok(PropStore { mmap, num_records: len / PROP_BYTES })
    }

    /// Number of property records in the store.
    #[inline]
    pub fn len(&self) -> usize {
        self.num_records
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.num_records == 0
    }

    /// Fetch the molecular properties for `doc_id`, or `None` if out of range.
    ///
    /// Bytes are copied out via `from_le_bytes`, so there are no alignment
    /// requirements on the underlying mmap.
    #[inline]
    pub fn get(&self, doc_id: u32) -> Option<MolecularProperties> {
        let idx = doc_id as usize;
        if idx >= self.num_records {
            return None;
        }
        let start = idx * PROP_BYTES;
        let bytes = &self.mmap[start..start + PROP_BYTES];
        let field = |off: usize| f32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
        Some(MolecularProperties {
            mw: field(0),
            logp: field(4),
            rot_bonds: field(8) as u32,
            heavy_atoms: field(12) as u32,
            ring_count: field(16) as u32,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::properties::compute_properties;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_store(props: &[MolecularProperties]) -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        for p in props {
            let mut buf = [0u8; PROP_BYTES];
            buf[0..4].copy_from_slice(&p.mw.to_le_bytes());
            buf[4..8].copy_from_slice(&p.logp.to_le_bytes());
            buf[8..12].copy_from_slice(&(p.rot_bonds as f32).to_le_bytes());
            buf[12..16].copy_from_slice(&(p.heavy_atoms as f32).to_le_bytes());
            buf[16..20].copy_from_slice(&(p.ring_count as f32).to_le_bytes());
            tmp.write_all(&buf).unwrap();
        }
        tmp.flush().unwrap();
        tmp
    }

    #[test]
    fn test_roundtrip() {
        let props = vec![
            compute_properties("CCO"),
            compute_properties("c1ccccc1"),
            compute_properties("CC(=O)Oc1ccccc1C(=O)O"),
        ];
        let tmp = write_store(&props);
        let store = PropStore::open(tmp.path()).unwrap();

        assert_eq!(store.len(), 3);
        for (i, p) in props.iter().enumerate() {
            let got = store.get(i as u32).unwrap();
            assert!((got.mw - p.mw).abs() < 1e-3);
            assert!((got.logp - p.logp).abs() < 1e-3);
            assert_eq!(got.rot_bonds, p.rot_bonds);
            assert_eq!(got.heavy_atoms, p.heavy_atoms);
            assert_eq!(got.ring_count, p.ring_count);
        }
        assert!(store.get(3).is_none());
    }
}
