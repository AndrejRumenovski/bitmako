//! Inverted index for 1024-bit Morgan fingerprints.

pub mod block_max;
pub mod builder;
pub mod posting_list;

use std::path::Path;

use memmap2::Mmap;
use tracing::info;

use crate::error::{BitMakoError, Result};
use crate::index::posting_list::PostingList;

const INDEX_MAGIC: &[u8; 8] = b"BITMAKO1";

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

        let _version = read_u32(8);
        let num_compounds = read_u32(12);
        let num_bits = read_u32(16);

        // Offset table: num_bits × u64, starting at byte 20.
        let offsets_start = 20usize;
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
        let start = self.postings_start + self.offsets[bit] as usize;
        let end = if bit + 1 < self.offsets.len() {
            self.postings_start + self.offsets[bit + 1] as usize
        } else {
            self.mmap.len()
        };
        PostingList::deserialize(&self.mmap[start..end]).map_err(BitMakoError::Io)
    }
}
