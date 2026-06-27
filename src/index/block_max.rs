//! Block-Max metadata for the BMW WAND engine.
//!
//! For each posting list (fingerprint bit), stores the maximum Tanimoto score
//! contribution possible from any document in each block.
//!
//! The maximum Tanimoto upper bound for a query term (bit b) in a block is:
//!
//!   UB(b, block) = 1 / max(popcount_Q, min_popcount_in_block)
//!
//! Where popcount_Q is the popcount of the query fingerprint, and
//! min_popcount_in_block is the minimum compound popcount in the block
//! (any compound with fewer bits set can have higher Tanimoto for this query bit).
//!
//! We pre-compute and store max_pop (the maximum compound popcount per block)
//! and derive upper bounds at query time.

use crate::index::posting_list::{PostingList, BLOCK_SIZE};

/// Per-bit, per-block upper bound on Tanimoto contribution.
///
/// `bounds[bit][block]` = max Tanimoto this bit can add from any doc in that block.
pub struct BlockMaxTable {
    /// Number of fingerprint bits (1024)
    pub num_bits: usize,
    /// For each bit, the max compound popcount per block (u8 = max 255 bits)
    /// Layout: [bit0_block0, bit0_block1, ..., bit1_block0, ...]
    pub max_pops: Vec<Vec<u8>>,
}

impl BlockMaxTable {
    pub fn new(num_bits: usize) -> Self {
        BlockMaxTable {
            num_bits,
            max_pops: vec![Vec::new(); num_bits],
        }
    }

    /// Set block-level max popcount for (bit, block)
    #[inline]
    pub fn set(&mut self, bit: usize, block: usize, max_pop: u8) {
        let pops = &mut self.max_pops[bit];
        if block >= pops.len() {
            pops.resize(block + 1, 0);
        }
        pops[block] = max_pop;
    }

    /// Compute the upper-bound Tanimoto score contribution from bit `bit`
    /// for documents in `block`, given that the query has `query_pop` bits set.
    ///
    /// Tanimoto(A, B) = |A ∩ B| / |A ∪ B| = |A ∩ B| / (|A| + |B| - |A ∩ B|)
    ///
    /// Upper bound when bit b is set in both query and document:
    ///   Each shared bit contributes 1 to the intersection.
    ///   The maximum total Tanimoto is bounded by:
    ///   intersection ≤ query_pop (can't share more bits than query has)
    ///
    /// For the WAND upper bound on a single-bit contribution:
    ///   UB ≈ 1 / (query_pop + doc_pop - intersection)
    ///   simplified to: 1 / max(query_pop, doc_min_pop_in_block)
    #[inline]
    pub fn upper_bound_tanimoto(&self, bit: usize, block: usize, query_pop: u32) -> f32 {
        let max_pop = self.max_pops[bit].get(block).copied().unwrap_or(0) as u32;
        if max_pop == 0 || query_pop == 0 {
            return 0.0;
        }
        // Upper bound: if the doc has max_pop bits set and query has query_pop bits set,
        // the intersection can be at most min(max_pop, query_pop).
        let max_intersection = max_pop.min(query_pop);
        let min_union = max_pop.max(query_pop); // union >= max(|A|, |B|)
        max_intersection as f32 / min_union as f32
    }

    /// Cumulative upper bound from all active bits in the query that have postings
    /// overlapping `block`. Used in BMW pivot selection.
    pub fn cumulative_upper_bound(
        &self,
        active_bits: &[usize],
        block: usize,
        query_pop: u32,
    ) -> f32 {
        active_bits
            .iter()
            .map(|&bit| self.upper_bound_tanimoto(bit, block, query_pop))
            .fold(0.0f32, |acc, x| acc + x)
            .min(1.0) // Tanimoto is always ≤ 1
    }
}

/// Build the BlockMaxTable from pre-built posting lists + per-compound popcounts.
///
/// `compound_pops[doc_id]` = popcount of that compound's fingerprint.
pub fn build_block_max_table(
    lists: &[PostingList],
    compound_pops: &[u8],
) -> BlockMaxTable {
    let mut table = BlockMaxTable::new(lists.len());

    for (bit, pl) in lists.iter().enumerate() {
        let num_blocks = pl.num_blocks();
        for block_idx in 0..num_blocks {
            let start = PostingList::block_start(block_idx);
            let end = (start + BLOCK_SIZE).min(pl.len());
            let max_pop = pl.doc_ids[start..end]
                .iter()
                .map(|&doc_id| compound_pops.get(doc_id as usize).copied().unwrap_or(0))
                .max()
                .unwrap_or(0);
            table.set(bit, block_idx, max_pop);
        }
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_upper_bound_tanimoto() {
        let mut table = BlockMaxTable::new(1);
        table.set(0, 0, 50); // block 0 of bit 0 has max_pop=50

        let ub = table.upper_bound_tanimoto(0, 0, 50);
        // max_intersection = min(50, 50) = 50, min_union = max(50,50) = 50
        // UB = 50/50 = 1.0
        assert!((ub - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_upper_bound_sparse_doc() {
        let mut table = BlockMaxTable::new(1);
        table.set(0, 0, 20); // doc in block has only 20 bits set
        // query has 100 bits set
        let ub = table.upper_bound_tanimoto(0, 0, 100);
        // max_intersection = min(20, 100) = 20, min_union = max(20,100) = 100
        // UB = 20/100 = 0.2
        assert!((ub - 0.2).abs() < 1e-6);
    }
}
