//! Block-Max WAND (BMW) similarity search engine.
//!
//! BMW implements dynamic pruning over inverted posting lists to find top-k
//! compounds by Tanimoto similarity without full-scan evaluation.
//!
//! Algorithm:
//!   1. Open a cursor over each posting list for active query bits.
//!   2. Advance in document-ID order, merging all cursors.
//!   3. For each candidate document:
//!      a. Apply block-level upper bound: if max Tanimoto possible in this
//!         block is below current threshold, jump to next block (block skip).
//!      b. Apply popcount upper bound: if popcount ratio cannot meet threshold,
//!         skip without computing full Tanimoto.
//!      c. Compute exact Tanimoto and update top-k heap.
//!   4. Dynamically raise threshold as heap fills, increasing skip rate.

use std::collections::BinaryHeap;
use std::cmp::Ordering;

use tracing::{debug};

use crate::error::Result;
use crate::etl::fingerprint::Fingerprint;
use crate::index::IndexReader;
use crate::index::posting_list::{PostingList, BLOCK_SIZE};
use crate::search::query::{SimilarityQuery, validate_query};
use crate::search::tanimoto::{tanimoto_upper_bound, tanimoto_with_threshold};

/// Scored candidate in the top-k min-heap.
/// BinaryHeap is a max-heap; we reverse comparison to make it a min-heap
/// so we efficiently evict the lowest-scoring element when heap overflows top_k.
#[derive(Debug, Clone)]
struct Candidate {
    score: ordered_float::OrderedFloat<f32>,
    doc_id: u32,
}

impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool { self.score == other.score }
}
impl Eq for Candidate {}
impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}
impl Ord for Candidate {
    // Reversed → min-heap: lowest score at the top for fast eviction
    fn cmp(&self, other: &Self) -> Ordering {
        other.score.cmp(&self.score)
    }
}

/// Iterator state over one posting list during search execution
struct ListCursor<'a> {
    #[allow(dead_code)]
    bit: usize,
    pl: &'a PostingList,
    pos: usize,
}

impl<'a> ListCursor<'a> {
    #[inline]
    fn current_doc(&self) -> Option<u32> {
        self.pl.doc_ids.get(self.pos).copied()
    }

    #[inline]
    fn is_exhausted(&self) -> bool {
        self.pos >= self.pl.len()
    }

    /// Advance to first position with doc_id >= target
    #[inline]
    fn advance_to(&mut self, target: u32) {
        self.pos = self.pl.advance_to(self.pos, target);
    }

    /// Max compound popcount in the current block (for block-level Tanimoto UB)
    #[inline]
    fn block_max_pop(&self) -> u8 {
        let block = self.pos / BLOCK_SIZE;
        self.pl.block_max(block)
    }

    /// Advance past the current block entirely
    #[inline]
    fn skip_block(&mut self) {
        let block = self.pos / BLOCK_SIZE;
        let next_block_start = (block + 1) * BLOCK_SIZE;
        self.pos = self.pl.advance_to(self.pos, next_block_start as u32);
    }
}

/// Block-Max WAND execution engine
pub struct BmwEngine<'a> {
    index: &'a IndexReader,
}

impl<'a> BmwEngine<'a> {
    pub fn new(index: &'a IndexReader) -> Self {
        BmwEngine { index }
    }

    /// Execute a BMW top-k Tanimoto similarity search.
    ///
    /// The `get_fingerprint` closure retrieves the full 1024-bit fingerprint for
    /// a given doc_id. In production this fetches from a memory-mapped flat file;
    /// in tests it can be a simple Vec lookup.
    ///
    /// Returns up to `query.top_k` results as `(doc_id, tanimoto_score)` pairs
    /// sorted by descending similarity.
    pub fn search(
        &self,
        query: &SimilarityQuery,
        get_fingerprint: impl Fn(u32) -> Option<Fingerprint>,
    ) -> Result<Vec<(u32, f32)>> {
        validate_query(query)?;

        let active_bits = query.active_bits();
        if active_bits.is_empty() {
            return Ok(Vec::new());
        }

        // Open one cursor per active bit posting list (skip empty lists)
        let mut cursors: Vec<ListCursor<'_>> = active_bits
            .iter()
            .filter_map(|&bit| {
                let pl = &self.index.posting_lists[bit];
                if pl.is_empty() {
                    None
                } else {
                    Some(ListCursor { bit, pl, pos: 0 })
                }
            })
            .collect();

        if cursors.is_empty() {
            return Ok(Vec::new());
        }

        let mut heap: BinaryHeap<Candidate> = BinaryHeap::with_capacity(query.top_k + 1);
        let mut threshold = query.threshold;
        let mut docs_evaluated: u64 = 0;
        let mut blocks_skipped: u64 = 0;
        let mut docs_skipped_by_popcount: u64 = 0;

        loop {
            // --- Step 1: Remove exhausted cursors ---
            cursors.retain(|c| !c.is_exhausted());
            if cursors.is_empty() {
                break;
            }

            // --- Step 2: Find minimum current doc_id (next candidate) ---
            let current_doc = cursors
                .iter()
                .filter_map(|c| c.current_doc())
                .min()
                .unwrap(); // safe: all cursors are non-exhausted

            // --- Step 3: Block-Max Upper Bound (BMW pruning) ---
            // For the cursor(s) at current_doc, check if ANY doc in their block
            // can possibly have Tanimoto >= threshold.
            let block_max_ub = self.compute_block_max_ub(query, &cursors, current_doc);

            if block_max_ub < threshold {
                // No document in this block can meet the threshold.
                // Advance all cursors that are in this block past it.
                blocks_skipped += 1;
                for cursor in cursors.iter_mut() {
                    if cursor.current_doc().map(|d| d <= current_doc).unwrap_or(false) {
                        cursor.skip_block();
                    }
                }
                continue;
            }

            // --- Step 4: Popcount upper bound (fast pre-filter) ---
            let candidate_pop = self
                .index
                .compound_pops
                .get(current_doc as usize)
                .copied()
                .unwrap_or(0) as u32;

            let pop_ub = tanimoto_upper_bound(query.query_pop, candidate_pop);

            if pop_ub < threshold {
                // Popcount ratio cannot meet threshold — skip this doc
                docs_skipped_by_popcount += 1;
                for cursor in cursors.iter_mut() {
                    if cursor.current_doc() == Some(current_doc) {
                        cursor.advance_to(current_doc + 1);
                    }
                }
                continue;
            }

            // --- Step 5: Exact Tanimoto evaluation ---
            if let Some(fp) = get_fingerprint(current_doc) {
                let (score, meets) = tanimoto_with_threshold(&query.query_fp, &fp, threshold);
                docs_evaluated += 1;

                if meets {
                    heap.push(Candidate {
                        score: ordered_float::OrderedFloat(score),
                        doc_id: current_doc,
                    });
                    if heap.len() > query.top_k {
                        heap.pop(); // evict lowest score
                    }
                    // Raise threshold to current worst kept score
                    if heap.len() == query.top_k {
                        if let Some(worst) = heap.peek() {
                            threshold = threshold.max(worst.score.0);
                        }
                    }
                }
            }

            // Advance all cursors that are at current_doc
            for cursor in cursors.iter_mut() {
                if cursor.current_doc() == Some(current_doc) {
                    cursor.advance_to(current_doc + 1);
                }
            }
        }

        debug!(
            "BMW stats: evaluated={} block_skipped={} pop_skipped={} results={}",
            docs_evaluated, blocks_skipped, docs_skipped_by_popcount, heap.len()
        );

        // `into_sorted_vec` on our min-heap (reversed Ord) yields elements in
        // ascending order of reversed-Ord == descending order of actual score.
        let results: Vec<(u32, f32)> = heap
            .into_sorted_vec()
            .into_iter()
            .map(|c| (c.doc_id, c.score.0))
            .collect();

        Ok(results)
    }

    /// Compute the block-level upper bound on Tanimoto.
    ///
    /// With only `max_block_pop` (highest compound popcount in the block):
    ///   - If query_pop ≤ max_block_pop: a doc with exactly query_pop bits might exist → UB = 1.0
    ///   - If query_pop > max_block_pop: all docs have fewer bits than the query → UB = max_pop / query_pop
    ///
    /// This is the tightest safe bound obtainable from max_block_pop alone.
    fn compute_block_max_ub(
        &self,
        query: &SimilarityQuery,
        cursors: &[ListCursor<'_>],
        current_doc: u32,
    ) -> f32 {
        let max_pop = cursors
            .iter()
            .filter(|c| c.current_doc().map(|d| d <= current_doc).unwrap_or(false))
            .map(|c| c.block_max_pop() as u32)
            .max()
            .unwrap_or(0);

        if max_pop == 0 || query.query_pop == 0 {
            return 0.0;
        }
        if query.query_pop <= max_pop {
            1.0
        } else {
            max_pop as f32 / query.query_pop as f32
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::fingerprint::compute_morgan_fp;
    use crate::index::builder::IndexBuilder;
    use crate::index::IndexReader;
    use crate::search::query::SimilarityQuery;
    use tempfile::NamedTempFile;

    fn build_index_from_smiles(smiles_list: &[&str]) -> (Vec<Fingerprint>, IndexReader, NamedTempFile) {
        let fps: Vec<Fingerprint> = smiles_list.iter().map(|s| compute_morgan_fp(s)).collect();
        let mut builder = IndexBuilder::new();
        for (i, fp) in fps.iter().enumerate() {
            builder.add_compound(i as u32, fp);
        }
        let tmp = NamedTempFile::new().unwrap();
        builder.write_index(tmp.path()).unwrap();
        let index = IndexReader::open(tmp.path()).unwrap();
        (fps, index, tmp)
    }

    #[test]
    fn test_exact_match_found() {
        let smiles = ["CCO", "c1ccccc1", "CNC(=O)c1ccccc1"];
        let (fps, index, _tmp) = build_index_from_smiles(&smiles);

        let query_fp = compute_morgan_fp("CCO");
        let query = SimilarityQuery::new(query_fp, 0.5, 5);
        let engine = BmwEngine::new(&index);

        let fps_c = fps.clone();
        let results = engine.search(&query, |doc_id| fps_c.get(doc_id as usize).copied()).unwrap();

        assert!(!results.is_empty());
        // CCO queried against itself should return score=1.0 as the top hit
        let (top_doc, top_score) = results[0];
        assert_eq!(top_doc, 0);
        assert!((top_score - 1.0).abs() < 1e-5, "Expected 1.0, got {}", top_score);
    }

    #[test]
    fn test_high_threshold_fewer_results() {
        let smiles = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "CN", "CC(=O)O"];
        let (fps, index, _tmp) = build_index_from_smiles(&smiles);
        let engine = BmwEngine::new(&index);
        let query_fp = compute_morgan_fp("CCO");

        let results_low = engine.search(
            &SimilarityQuery::new(query_fp, 0.3, 100),
            |doc_id| fps.get(doc_id as usize).copied()
        ).unwrap();

        let results_high = engine.search(
            &SimilarityQuery::new(query_fp, 0.9, 100),
            |doc_id| fps.get(doc_id as usize).copied()
        ).unwrap();

        // Higher threshold → fewer or equal results
        assert!(results_high.len() <= results_low.len());
    }

    #[test]
    fn test_results_sorted_descending() {
        let smiles = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "CN", "CC(=O)O"];
        let (fps, index, _tmp) = build_index_from_smiles(&smiles);
        let engine = BmwEngine::new(&index);
        let query_fp = compute_morgan_fp("CCO");
        let query = SimilarityQuery::new(query_fp, 0.0, 100);

        let results = engine.search(&query, |doc_id| fps.get(doc_id as usize).copied()).unwrap();

        // Results from into_sorted_vec on a reversed-Ord heap are descending by score
        for window in results.windows(2) {
            assert!(
                window[0].1 >= window[1].1,
                "Results not sorted descending: {} < {}",
                window[0].1, window[1].1
            );
        }
    }
}
