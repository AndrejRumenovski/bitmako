//! Inverted index for 1024-bit Morgan fingerprints.

pub mod block_max;
pub mod builder;
pub mod posting_list;
pub mod skip;

use std::path::Path;

use memmap2::Mmap;
use tracing::info;

use crate::error::{BitMakoError, Result};
use crate::index::posting_list::PostingList;

/// On-disk index file magic bytes, shared by every writer of the `*.bitmako`
/// format (`IndexBuilder::write_index` and the streaming multi-pass builder in
/// `main.rs`'s `cmd_build_index`) and validated here on read.
pub const INDEX_MAGIC: &[u8; 8] = b"BITMAKO1";
/// Current on-disk index format version, written by every builder as of this
/// commit: 8-byte (u64) posting-list offsets, supporting posting sections over
/// 4 GiB.
///
/// Note: an earlier build of `cmd_build_index` wrote version `1` into the header
/// by mistake while *already* using this same 8-byte-offset layout (a 4-byte-offset
/// v1 format was never actually implemented anywhere in this codebase). `IndexReader`
/// therefore accepts both `1` and `2` on read — see `SUPPORTED_INDEX_VERSIONS` — so
/// existing index files built before this fix keep working.
pub const INDEX_VERSION: u32 = 2;
/// Version values `IndexReader::open` accepts; see `INDEX_VERSION`.
const SUPPORTED_INDEX_VERSIONS: [u32; 2] = [1, 2];

/// Memory-mapped read-only view of a serialized inverted index.
///
/// The index is **not** fully decoded on open — at 1.4B compounds the decoded
/// posting lists would need ~190 GB of RAM. Instead the file is mmap'd and:
///   - `compound_pop(doc_id)` reads a single byte directly from the map,
///   - `decode_posting_list(bit)` decodes exactly one posting list on demand.
///
/// A query touches only the posting lists for its active bits, so peak memory
/// stays proportional to the query, not the corpus.
pub struct IndexReader {
    mmap: Mmap,
    pub num_compounds: u32,
    pub num_bits: u32,
    /// Byte offset of each posting list relative to `postings_start`.
    offsets: Vec<u64>,
    /// Byte offset where the compound popcount table begins.
    pops_start: usize,
    /// Byte offset where the posting list data section begins.
    postings_start: usize,
}

impl IndexReader {
    /// Open a memory-mapped index file built by IndexBuilder::write_index.
    pub fn open(path: &Path) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Self::from_mmap(mmap)
    }

    fn from_mmap(mmap: Mmap) -> Result<Self> {
        let data = &mmap[..];
        if data.len() < 20 || &data[0..8] != INDEX_MAGIC {
            return Err(BitMakoError::IndexBuild("invalid index magic bytes".into()));
        }

        let read_u32 = |off: usize| -> u32 {
            u32::from_le_bytes(data[off..off + 4].try_into().unwrap())
        };

        let version = read_u32(8);
        if !SUPPORTED_INDEX_VERSIONS.contains(&version) {
            return Err(BitMakoError::IndexBuild(format!(
                "unsupported index version {} (expected one of {:?})",
                version, SUPPORTED_INDEX_VERSIONS
            )));
        }
        let num_compounds = read_u32(12);
        let num_bits = read_u32(16);

        // Offset table: num_bits × u64, starting at byte 20.
        let offsets_start = 20usize;
        let offsets_end = offsets_start + num_bits as usize * 8;
        if offsets_end > data.len() {
            return Err(BitMakoError::IndexBuild(
                "index truncated: offset table exceeds file length".into(),
            ));
        }
        let mut offsets = vec![0u64; num_bits as usize];
        for (i, off) in offsets.iter_mut().enumerate() {
            let p = offsets_start + i * 8;
            *off = u64::from_le_bytes(data[p..p + 8].try_into().unwrap());
        }

        // Compound pops follow the offset table; postings follow pops + padding.
        let pops_start = offsets_start + num_bits as usize * 8;
        let pad = (4 - (num_compounds as usize % 4)) % 4;
        let postings_start = pops_start + num_compounds as usize + pad;

        if postings_start > data.len() {
            return Err(BitMakoError::IndexBuild(
                "index truncated: header exceeds file length".into(),
            ));
        }

        info!(
            "Opened index (mmap): {} compounds, {} bits, {} MB",
            num_compounds,
            num_bits,
            data.len() / (1024 * 1024)
        );

        Ok(IndexReader {
            mmap,
            num_compounds,
            num_bits,
            offsets,
            pops_start,
            postings_start,
        })
    }

    /// Popcount of the fingerprint for `doc_id` (0 if out of range).
    #[inline]
    pub fn compound_pop(&self, doc_id: u32) -> u8 {
        self.mmap
            .get(self.pops_start + doc_id as usize)
            .copied()
            .unwrap_or(0)
    }

    /// Decode the posting list for a single fingerprint bit.
    ///
    /// Only the bytes for this one list are read from the mmap and decoded,
    /// so callers should decode just the active bits of their query.
    pub fn decode_posting_list(&self, bit: usize) -> Result<PostingList> {
        PostingList::deserialize(self.posting_bytes(bit)).map_err(BitMakoError::Io)
    }

    /// Raw serialized bytes of the posting list for `bit` (the `u32 num_blocks`
    /// header followed by block data). Used by the streaming cursor and the skip
    /// index builder, which decode block-by-block rather than all at once.
    #[inline]
    pub fn posting_bytes(&self, bit: usize) -> &[u8] {
        let start = self.postings_start + self.offsets[bit] as usize;
        let end = if bit + 1 < self.offsets.len() {
            self.postings_start + self.offsets[bit + 1] as usize
        } else {
            self.mmap.len()
        };
        &self.mmap[start..end]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::fingerprint::compute_morgan_fp;
    use crate::index::builder::IndexBuilder;

    fn write_and_open(smiles: &[&str]) -> (IndexReader, tempfile::NamedTempFile) {
        let mut builder = IndexBuilder::new();
        for (i, s) in smiles.iter().enumerate() {
            builder.add_compound(i as u32, &compute_morgan_fp(s));
        }
        let tmp = tempfile::NamedTempFile::new().unwrap();
        builder.write_index(tmp.path()).unwrap();
        let reader = IndexReader::open(tmp.path()).unwrap();
        (reader, tmp)
    }

    #[test]
    fn open_reports_correct_header_fields() {
        let (reader, _tmp) = write_and_open(&["CCO", "c1ccccc1", "CNC(=O)c1ccccc1"]);
        assert_eq!(reader.num_compounds, 3);
        assert_eq!(reader.num_bits, crate::etl::fingerprint::FP_BITS as u32);
    }

    #[test]
    fn compound_pop_matches_actual_fingerprint_popcount() {
        let smiles = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"];
        let (reader, _tmp) = write_and_open(&smiles);
        for (i, s) in smiles.iter().enumerate() {
            let fp = compute_morgan_fp(s);
            let expected_pop = crate::etl::fingerprint::fp_popcount(&fp) as u8;
            assert_eq!(reader.compound_pop(i as u32), expected_pop);
        }
    }

    #[test]
    fn compound_pop_out_of_range_returns_zero() {
        let (reader, _tmp) = write_and_open(&["CCO"]);
        assert_eq!(reader.compound_pop(9999), 0);
    }

    #[test]
    fn decode_posting_list_contains_every_doc_with_that_bit_set() {
        let smiles = ["CCO", "CCCO", "c1ccccc1", "CN"];
        let (reader, _tmp) = write_and_open(&smiles);
        let fps: Vec<_> = smiles.iter().map(|s| compute_morgan_fp(s)).collect();

        // Pick a bit that's actually set in at least one fingerprint and verify
        // the decoded posting list is exactly the set of docs with that bit on.
        for bit in 0..crate::etl::fingerprint::FP_BITS {
            let word = bit / 64;
            let b = bit % 64;
            let expected: Vec<u32> = fps
                .iter()
                .enumerate()
                .filter(|(_, fp)| (fp[word] >> b) & 1 == 1)
                .map(|(i, _)| i as u32)
                .collect();
            if expected.is_empty() {
                continue;
            }
            let decoded = reader.decode_posting_list(bit).unwrap();
            assert_eq!(decoded.doc_ids, expected, "mismatch for bit {bit}");
        }
    }

    #[test]
    fn open_rejects_bad_magic_bytes() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"NOTAMAGIC00000000000").unwrap();
        assert!(IndexReader::open(tmp.path()).is_err());
    }

    #[test]
    fn open_rejects_unsupported_version() {
        let mut data = Vec::new();
        data.extend_from_slice(INDEX_MAGIC);
        data.extend_from_slice(&99u32.to_le_bytes()); // version 99 unsupported
        data.extend_from_slice(&0u32.to_le_bytes()); // num_compounds
        data.extend_from_slice(&0u32.to_le_bytes()); // num_bits
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), &data).unwrap();
        assert!(IndexReader::open(tmp.path()).is_err());
    }

    #[test]
    fn open_accepts_legacy_version_1_header() {
        // Historical bug: cmd_build_index once wrote version byte 1 while
        // already using the 8-byte-offset v2 layout. IndexReader must keep
        // accepting real files built with that header — see INDEX_VERSION's
        // doc comment. Build a normal index then patch the version field.
        let mut builder = IndexBuilder::new();
        builder.add_compound(0, &compute_morgan_fp("CCO"));
        let tmp = tempfile::NamedTempFile::new().unwrap();
        builder.write_index(tmp.path()).unwrap();

        let mut bytes = std::fs::read(tmp.path()).unwrap();
        bytes[8..12].copy_from_slice(&1u32.to_le_bytes());
        std::fs::write(tmp.path(), &bytes).unwrap();

        let reader = IndexReader::open(tmp.path()).unwrap();
        assert_eq!(reader.num_compounds, 1);
    }

    #[test]
    fn open_rejects_truncated_header() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        // Magic bytes only, well under the 20-byte minimum header.
        std::fs::write(tmp.path(), INDEX_MAGIC).unwrap();
        assert!(IndexReader::open(tmp.path()).is_err());
    }

    #[test]
    fn open_rejects_offset_table_exceeding_file_length() {
        // A header claiming a huge num_bits (so its offsets table alone would
        // exceed the tiny file we actually write) must error, not panic on
        // out-of-bounds slicing.
        let mut data = Vec::new();
        data.extend_from_slice(INDEX_MAGIC);
        data.extend_from_slice(&INDEX_VERSION.to_le_bytes());
        data.extend_from_slice(&1u32.to_le_bytes()); // num_compounds = 1
        data.extend_from_slice(&1024u32.to_le_bytes()); // num_bits = 1024
        // Deliberately omit the offsets table / pops / postings sections.
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), &data).unwrap();
        assert!(IndexReader::open(tmp.path()).is_err());
    }
}
