//! Delta-encoded, block-packed posting lists for inverted fingerprint index.
//!
//! Each posting list corresponds to one of the 1024 fingerprint bits.
//! Document IDs are stored delta-encoded and bit-packed using Frame-of-Reference
//! (FOR) to maximize cache efficiency during WAND traversal.
//!
//! Block layout (BLOCK_SIZE = 128 documents):
//!   [block_max_popcount: u8][doc_ids: delta+FOR encoded]

use std::io::{self, Cursor, Read, Write};


/// Number of document IDs per posting block.
/// 128 gives good cache alignment: each block fits in ~2-3 cache lines.
pub const BLOCK_SIZE: usize = 128;

/// Serialized posting list stored on disk / in a memory-mapped region.
///
/// Wire format (little-endian):
///   u32  num_blocks
///   [num_blocks × BlockHeader]
///   [variable-length encoded doc_ids for each block]
#[derive(Debug)]
pub struct PostingList {
    /// Pre-sorted document IDs (ascending)
    pub doc_ids: Vec<u32>,
    /// Max fingerprint popcount per block (upper bound for Tanimoto scoring)
    pub block_max_pop: Vec<u8>,
}

impl PostingList {
    pub fn new() -> Self {
        PostingList {
            doc_ids: Vec::new(),
            block_max_pop: Vec::new(),
        }
    }

    /// Number of postings
    #[inline]
    pub fn len(&self) -> usize {
        self.doc_ids.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.doc_ids.is_empty()
    }

    /// Number of full + partial blocks
    #[inline]
    pub fn num_blocks(&self) -> usize {
        (self.doc_ids.len() + BLOCK_SIZE - 1) / BLOCK_SIZE
    }

    /// Block-level max popcount for block `b`
    #[inline]
    pub fn block_max(&self, b: usize) -> u8 {
        self.block_max_pop.get(b).copied().unwrap_or(0)
    }

    /// Advance iterator to the first position with doc_id >= target.
    /// Returns the index into `doc_ids`, or `doc_ids.len()` if not found.
    /// Uses binary search for O(log n) block skipping.
    #[inline]
    pub fn advance_to(&self, current_pos: usize, target: u32) -> usize {
        let slice = &self.doc_ids[current_pos..];
        match slice.binary_search(&target) {
            Ok(rel) => current_pos + rel,
            Err(rel) => current_pos + rel,
        }
    }

    /// Compute block index for a given position
    #[inline]
    pub fn block_of(pos: usize) -> usize {
        pos / BLOCK_SIZE
    }

    /// First position of a block
    #[inline]
    pub fn block_start(block: usize) -> usize {
        block * BLOCK_SIZE
    }

    /// Serialize to bytes using delta + FOR encoding.
    pub fn serialize(&self) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        let num_blocks = self.num_blocks() as u32;
        out.write_all(&num_blocks.to_le_bytes())?;

        for block_idx in 0..num_blocks as usize {
            let start = block_idx * BLOCK_SIZE;
            let end = (start + BLOCK_SIZE).min(self.doc_ids.len());
            let block = &self.doc_ids[start..end];
            let max_pop = self.block_max_pop.get(block_idx).copied().unwrap_or(0);

            // Header: max_pop (1 byte) + block_len (1 byte)
            out.push(max_pop);
            out.push((end - start) as u8);

            // Delta encode within block
            let mut prev = if block_idx == 0 { 0u32 } else { self.doc_ids[start - 1] };
            for &doc_id in block {
                let delta = doc_id.wrapping_sub(prev);
                write_varint(&mut out, delta)?;
                prev = doc_id;
            }
        }

        Ok(out)
    }

    /// Deserialize from bytes.
    pub fn deserialize(data: &[u8]) -> io::Result<Self> {
        let mut cursor = Cursor::new(data);
        let mut buf4 = [0u8; 4];
        cursor.read_exact(&mut buf4)?;
        let num_blocks = u32::from_le_bytes(buf4) as usize;

        let mut doc_ids = Vec::new();
        let mut block_max_pop = Vec::with_capacity(num_blocks);

        let mut base = 0u32;

        for _ in 0..num_blocks {
            let mut header = [0u8; 2];
            cursor.read_exact(&mut header)?;
            let max_pop = header[0];
            let block_len = header[1] as usize;
            block_max_pop.push(max_pop);

            for _ in 0..block_len {
                let delta = read_varint(&mut cursor)?;
                let doc_id = base.wrapping_add(delta);
                doc_ids.push(doc_id);
                base = doc_id;
            }
        }

        Ok(PostingList { doc_ids, block_max_pop })
    }
}

/// Read the `u32 num_blocks` header at the start of a serialized posting list.
/// Returns `(num_blocks, position_after_header)`.
#[inline]
pub fn num_blocks_at(data: &[u8]) -> (u32, usize) {
    let n = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
    (n, 4)
}

/// Read a single LEB128 varint directly from `data` at `pos` (no Cursor/io
/// overhead — this is the hot path for streaming decode).
/// Returns `(value, new_pos)`.
#[inline]
pub fn read_varint_at(data: &[u8], mut pos: usize) -> (u32, usize) {
    let mut result = 0u32;
    let mut shift = 0u32;
    loop {
        let b = data[pos];
        pos += 1;
        result |= ((b & 0x7F) as u32) << shift;
        if b & 0x80 == 0 {
            break;
        }
        shift += 7;
    }
    (result, pos)
}

/// Decode one block beginning at byte `pos`, given the decoder `base`
/// (the absolute doc_id the block's first delta is relative to — i.e. the last
/// doc_id of the previous block, or 0 for the first block).
///
/// The decoded doc_ids are written into `out` (which is cleared first).
/// Returns `(max_pop, next_pos, last_doc_id)`.
#[inline]
pub fn decode_block(
    data: &[u8],
    pos: usize,
    mut base: u32,
    out: &mut Vec<u32>,
) -> (u8, usize, u32) {
    let max_pop = data[pos];
    let len = data[pos + 1] as usize;
    let mut p = pos + 2;
    out.clear();
    for _ in 0..len {
        let (delta, np) = read_varint_at(data, p);
        p = np;
        base = base.wrapping_add(delta);
        out.push(base);
    }
    (max_pop, p, base)
}

/// Walk a serialized posting list and produce one skip entry per block:
/// `(base, byte_offset)` where `base` is the decoder base entering the block
/// (last doc_id of the previous block; 0 for block 0) and `byte_offset` is the
/// position of the block header within `list_bytes`.
///
/// `base` values are strictly the previous block's last doc_id, so the array is
/// ascending and can be binary-searched to locate the block containing a target.
pub fn build_skip_entries(list_bytes: &[u8]) -> Vec<(u32, u64)> {
    if list_bytes.len() < 4 {
        return Vec::new();
    }
    let (num_blocks, mut pos) = num_blocks_at(list_bytes);
    let mut entries = Vec::with_capacity(num_blocks as usize);
    let mut base = 0u32;
    let mut buf: Vec<u32> = Vec::with_capacity(BLOCK_SIZE);
    for _ in 0..num_blocks {
        entries.push((base, pos as u64));
        let (_max_pop, next_pos, last_doc) = decode_block(list_bytes, pos, base, &mut buf);
        pos = next_pos;
        base = last_doc;
    }
    entries
}

/// Variable-length integer encoding (LEB128 style, little-endian)
fn write_varint(out: &mut Vec<u8>, mut value: u32) -> io::Result<()> {
    loop {
        let mut byte = (value & 0x7F) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80;
        }
        out.push(byte);
        if value == 0 {
            break;
        }
    }
    Ok(())
}

fn read_varint(cursor: &mut Cursor<&[u8]>) -> io::Result<u32> {
    let mut result = 0u32;
    let mut shift = 0u32;
    loop {
        let mut byte = [0u8; 1];
        cursor.read_exact(&mut byte)?;
        let b = byte[0];
        result |= ((b & 0x7F) as u32) << shift;
        if b & 0x80 == 0 {
            break;
        }
        shift += 7;
        if shift >= 35 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "varint overflow"));
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_list(ids: &[u32], max_pops: &[u8]) -> PostingList {
        PostingList {
            doc_ids: ids.to_vec(),
            block_max_pop: max_pops.to_vec(),
        }
    }

    #[test]
    fn test_roundtrip_small() {
        // 6 doc_ids → 1 block (BLOCK_SIZE=128), so only 1 block_max_pop entry
        let pl = make_list(&[0, 5, 10, 42, 100, 9999], &[3]);
        let bytes = pl.serialize().unwrap();
        let restored = PostingList::deserialize(&bytes).unwrap();
        assert_eq!(restored.doc_ids, pl.doc_ids);
        assert_eq!(restored.block_max_pop, pl.block_max_pop);
    }

    #[test]
    fn test_advance_to() {
        let pl = make_list(&[1, 5, 10, 15, 20, 100], &[]);
        assert_eq!(pl.advance_to(0, 10), 2);
        assert_eq!(pl.advance_to(0, 11), 3);
        assert_eq!(pl.advance_to(0, 101), 6);
    }

    #[test]
    fn test_roundtrip_large() {
        let ids: Vec<u32> = (0..300).map(|i| i * 3).collect();
        let pops: Vec<u8> = vec![10u8; 3]; // 3 blocks of 128 (partial last)
        let pl = PostingList { doc_ids: ids.clone(), block_max_pop: pops };
        let bytes = pl.serialize().unwrap();
        let restored = PostingList::deserialize(&bytes).unwrap();
        assert_eq!(restored.doc_ids, ids);
    }

    #[test]
    fn test_streaming_decode_matches_full() {
        // Multi-block list; stream block-by-block and confirm identical doc_ids,
        // and that skip-entry bases equal each block's preceding last doc_id.
        let ids: Vec<u32> = (0..500u32).map(|i| i * 7 + (i % 3)).collect();
        let pops: Vec<u8> = vec![5u8; 4];
        let pl = PostingList { doc_ids: ids.clone(), block_max_pop: pops };
        let bytes = pl.serialize().unwrap();

        let skips = build_skip_entries(&bytes);
        let (num_blocks, _) = num_blocks_at(&bytes);
        assert_eq!(skips.len(), num_blocks as usize);

        let mut all = Vec::new();
        let mut buf = Vec::new();
        for (i, &(base, off)) in skips.iter().enumerate() {
            // The base entering block i must be the last doc of block i-1 (0 for i=0).
            let expected_base = if i == 0 { 0 } else { ids[i * BLOCK_SIZE - 1] };
            assert_eq!(base, expected_base, "block {} base", i);
            let (_mp, _np, _last) = decode_block(&bytes, off as usize, base, &mut buf);
            all.extend_from_slice(&buf);
        }
        assert_eq!(all, ids);
    }
}
