//! Skip index ("sidecar") for streaming posting-list traversal.
//!
//! The inverted index stores doc_ids delta+varint encoded in 128-doc blocks with
//! deltas continuous *across* blocks, so a list can't be entered mid-stream
//! without knowing the running base. The skip index records, for every block of
//! every bit, the decoder base entering the block and the byte offset of the
//! block header. With it, a cursor can binary-search to the block containing a
//! target doc_id and decode only that block — turning `advance_to` from O(list)
//! into O(log blocks + one block), which keeps common-fragment queries from
//! decoding (or OOM-ing on) hundreds of millions of postings.
//!
//! On-disk layout (little-endian):
//! ```text
//! [8B magic "BMSKIP01"][4B version][4B num_bits]
//! [num_bits × { num_blocks: u32, data_off: u64 }]   // directory
//! [per bit: num_blocks × { base: u32, byte_offset: u64 }]
//! ```
//! `byte_offset` is relative to the start of that bit's `posting_bytes` slice.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use memmap2::Mmap;

use crate::error::{BitMakoError, Result};
use crate::index::posting_list::{build_skip_entries, num_blocks_at};
use crate::index::IndexReader;

const SKIP_MAGIC: &[u8; 8] = b"BMSKIP01";
const HEADER_LEN: usize = 8 + 4 + 4;
const DIR_ENTRY_LEN: usize = 12; // u32 num_blocks + u64 data_off
const SKIP_ENTRY_LEN: usize = 12; // u32 base + u64 byte_offset

enum Backing {
    Owned(Vec<u8>),
    Mapped(Mmap),
}

impl Backing {
    #[inline]
    fn bytes(&self) -> &[u8] {
        match self {
            Backing::Owned(v) => v,
            Backing::Mapped(m) => m,
        }
    }
}

/// Read-only skip index over a sidecar file (mmap) or in-memory bytes (tests).
pub struct SkipIndex {
    backing: Backing,
    num_bits: u32,
}

impl SkipIndex {
    /// Serialize a skip index for `index` into any writer.
    fn serialize<W: Write>(index: &IndexReader, w: &mut W) -> io::Result<()> {
        let num_bits = index.num_bits as usize;

        // Pass A: number of blocks per bit (cheap — just the list's u32 header).
        let mut nblocks = vec![0u32; num_bits];
        for (bit, nb) in nblocks.iter_mut().enumerate() {
            let lb = index.posting_bytes(bit);
            *nb = if lb.len() >= 4 { num_blocks_at(lb).0 } else { 0 };
        }

        // Header.
        w.write_all(SKIP_MAGIC)?;
        w.write_all(&1u32.to_le_bytes())?;
        w.write_all(&(num_bits as u32).to_le_bytes())?;

        // Directory: (num_blocks, data_off) per bit.
        let mut data_off = (HEADER_LEN + num_bits * DIR_ENTRY_LEN) as u64;
        for &nb in &nblocks {
            w.write_all(&nb.to_le_bytes())?;
            w.write_all(&data_off.to_le_bytes())?;
            data_off += nb as u64 * SKIP_ENTRY_LEN as u64;
        }

        // Entries: decode each list once, emit (base, byte_offset) per block.
        for (bit, &nb) in nblocks.iter().enumerate() {
            let lb = index.posting_bytes(bit);
            let entries = build_skip_entries(lb);
            debug_assert_eq!(entries.len(), nb as usize);
            for (base, off) in entries {
                w.write_all(&base.to_le_bytes())?;
                w.write_all(&off.to_le_bytes())?;
            }
        }
        Ok(())
    }

    /// Build the skip index for `index` and write it to `path`.
    pub fn build_to_file(index: &IndexReader, path: &Path) -> Result<()> {
        let file = File::create(path)?;
        let mut w = BufWriter::with_capacity(8 * 1024 * 1024, file);
        Self::serialize(index, &mut w).map_err(BitMakoError::Io)?;
        w.flush().map_err(BitMakoError::Io)?;
        Ok(())
    }

    /// Build the skip index entirely in memory (used by tests).
    pub fn build_in_memory(index: &IndexReader) -> Result<Self> {
        let mut buf = Vec::new();
        Self::serialize(index, &mut buf).map_err(BitMakoError::Io)?;
        Self::from_backing(Backing::Owned(buf))
    }

    /// Open a skip index file via mmap.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Self::from_backing(Backing::Mapped(mmap))
    }

    fn from_backing(backing: Backing) -> Result<Self> {
        let data = backing.bytes();
        if data.len() < HEADER_LEN || &data[0..8] != SKIP_MAGIC {
            return Err(BitMakoError::IndexBuild("invalid skip index magic".into()));
        }
        let num_bits = u32::from_le_bytes(data[12..16].try_into().unwrap());
        if data.len() < HEADER_LEN + num_bits as usize * DIR_ENTRY_LEN {
            return Err(BitMakoError::IndexBuild("skip index truncated directory".into()));
        }
        Ok(SkipIndex { backing, num_bits })
    }

    #[inline]
    pub fn num_bits(&self) -> u32 {
        self.num_bits
    }

    /// Skip entries for one bit's posting list.
    #[inline]
    pub fn entries(&self, bit: usize) -> SkipSlice<'_> {
        let data = self.backing.bytes();
        let dir = HEADER_LEN + bit * DIR_ENTRY_LEN;
        let num_blocks = u32::from_le_bytes(data[dir..dir + 4].try_into().unwrap()) as usize;
        let data_off = u64::from_le_bytes(data[dir + 4..dir + 12].try_into().unwrap()) as usize;
        SkipSlice {
            data: &data[data_off..data_off + num_blocks * SKIP_ENTRY_LEN],
            num_blocks,
        }
    }
}

/// A borrowed view of one bit's skip entries, supporting block lookup by doc_id.
pub struct SkipSlice<'a> {
    data: &'a [u8],
    num_blocks: usize,
}

impl<'a> SkipSlice<'a> {
    #[inline]
    pub fn num_blocks(&self) -> usize {
        self.num_blocks
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.num_blocks == 0
    }

    /// Decoder base entering block `i` (last doc_id of block `i-1`; 0 for block 0).
    #[inline]
    pub fn base(&self, i: usize) -> u32 {
        let off = i * SKIP_ENTRY_LEN;
        u32::from_le_bytes(self.data[off..off + 4].try_into().unwrap())
    }

    /// Byte offset of block `i`'s header within the bit's posting bytes.
    #[inline]
    pub fn byte_offset(&self, i: usize) -> usize {
        let off = i * SKIP_ENTRY_LEN + 4;
        u64::from_le_bytes(self.data[off..off + 8].try_into().unwrap()) as usize
    }

    /// Index of the block that may contain `target`: the largest block whose
    /// entering base is < target (clamped to block 0). Because base values are
    /// the previous block's last doc_id, the chosen block's doc range
    /// `(base, last_doc]` is guaranteed to cover `target` unless `target`
    /// exceeds the whole list (handled by the caller when scanning finds nothing).
    #[inline]
    pub fn block_for(&self, target: u32) -> usize {
        // Binary search for the count of blocks with base < target.
        let mut lo = 0usize;
        let mut hi = self.num_blocks;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.base(mid) < target {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        // `lo` is the count with base < target; the containing block is lo-1.
        lo.saturating_sub(1)
    }
}
