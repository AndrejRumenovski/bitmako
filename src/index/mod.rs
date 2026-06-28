//! Inverted index for 1024-bit Morgan fingerprints.

pub mod block_max;
pub mod builder;
pub mod posting_list;

use std::io::{Cursor, Read};
use std::path::Path;

use memmap2::Mmap;
use tracing::info;

use crate::error::{BitMakoError, Result};
use crate::index::posting_list::PostingList;

const INDEX_MAGIC: &[u8; 8] = b"BITMAKO1";

/// Memory-mapped read-only view of a serialized inverted index.
pub struct IndexReader {
    pub num_compounds: u32,
    pub num_bits: u32,
    pub compound_pops: Vec<u8>,
    pub posting_lists: Vec<PostingList>,
}

impl IndexReader {
    /// Open a memory-mapped index file built by IndexBuilder::write_index.
    pub fn open(path: &Path) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Self::from_bytes(&mmap)
    }

    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut cursor = Cursor::new(data);

        // Magic
        let mut magic = [0u8; 8];
        cursor.read_exact(&mut magic).map_err(BitMakoError::Io)?;
        if &magic != INDEX_MAGIC {
            return Err(BitMakoError::IndexBuild("invalid index magic bytes".into()));
        }

        let mut buf4 = [0u8; 4];

        cursor.read_exact(&mut buf4).map_err(BitMakoError::Io)?;
        let _version = u32::from_le_bytes(buf4);

        cursor.read_exact(&mut buf4).map_err(BitMakoError::Io)?;
        let num_compounds = u32::from_le_bytes(buf4);

        cursor.read_exact(&mut buf4).map_err(BitMakoError::Io)?;
        let num_bits = u32::from_le_bytes(buf4);

        // Read offsets (u64, 8 bytes each)
        let mut offsets = vec![0u64; num_bits as usize];
        for off in offsets.iter_mut() {
            let mut buf8 = [0u8; 8];
            cursor.read_exact(&mut buf8).map_err(BitMakoError::Io)?;
            *off = u64::from_le_bytes(buf8);
        }

        // Read compound pops
        let mut compound_pops = vec![0u8; num_compounds as usize];
        cursor.read_exact(&mut compound_pops).map_err(BitMakoError::Io)?;
        // Skip padding
        let pad = (4 - (num_compounds as usize % 4)) % 4;
        for _ in 0..pad {
            let mut _b = [0u8; 1];
            cursor.read_exact(&mut _b).map_err(BitMakoError::Io)?;
        }

        // Read posting lists from remaining bytes
        let postings_start = cursor.position() as usize;
        let postings_data = &data[postings_start..];

        let mut posting_lists = Vec::with_capacity(num_bits as usize);
        for i in 0..num_bits as usize {
            let start = offsets[i] as usize;
            let end = if i + 1 < offsets.len() {
                offsets[i + 1] as usize
            } else {
                postings_data.len()
            };
            let pl = PostingList::deserialize(&postings_data[start..end])
                .map_err(BitMakoError::Io)?;
            posting_lists.push(pl);
        }


        info!(
            "Loaded index: {} compounds, {} bits",
            num_compounds, num_bits
        );

        Ok(IndexReader {
            num_compounds,
            num_bits,
            compound_pops,
            posting_lists,
        })
    }
}
