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
use crate::index::posting_list::PostingList;
use crate::search::query::{SimilarityQuery, validate_query};
use crate::search::tanimoto::tanimoto_with_threshold;

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

        // Decode only the posting lists for the query's active bits. At corpus
        // scale this is the difference between a few GB and ~190 GB of RAM.
        let decoded: Vec<(usize, PostingList)> = active_bits
            .iter()
            .map(|&bit| Ok::<_, crate::error::BitMakoError>((bit, self.index.decode_posting_list(bit)?)))
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .filter(|(_, pl)| !pl.is_empty())
            .collect();

        // Open one cursor per non-empty active-bit posting list.
        let mut cursors: Vec<ListCursor<'_>> = decoded
            .iter()
            .map(|(bit, pl)| ListCursor { bit: *bit, pl, pos: 0 })
            .collect();

        if cursors.is_empty() {
            return Ok(Vec::new());
        }

        let mut heap: BinaryHeap<Candidate> = BinaryHeap::with_capacity(query.top_k + 1);
        let mut threshold = query.threshold;
        let p = query.query_pop;
        let mut docs_evaluated: u64 = 0;
        let mut pivots_skipped: u64 = 0;

        // Minimum number of query bits a candidate must share to possibly reach
        // `t`. A doc sharing `c` query bits has popcount K ≥ c, so its Tanimoto is
        // c/(P+K−c) ≤ c/P; reaching `t` therefore requires c ≥ t·P. This is the
        // WAND pivot threshold — only docs appearing in ≥ θ posting lists can win.
        let min_shared = |t: f32| -> usize {
            if t <= 0.0 {
                1
            } else {
                ((t * p as f32).ceil() as usize).clamp(1, p as usize)
            }
        };

        loop {
            cursors.retain(|c| !c.is_exhausted());
            let theta = min_shared(threshold);
            if theta > cursors.len() {
                break; // not enough remaining lists to ever reach θ shared bits
            }

            // Sort live cursors by current doc; the θ-th smallest is the pivot.
            cursors.sort_unstable_by_key(|c| c.current_doc().unwrap());
            let pivot_doc = cursors[theta - 1].current_doc().unwrap();

            if cursors[0].current_doc().unwrap() == pivot_doc {
                // Fully aligned: every cursor at `pivot_doc` is a shared query bit,
                // so the count gives the exact intersection size c.
                let c = cursors
                    .iter()
                    .filter(|cur| cur.current_doc() == Some(pivot_doc))
                    .count() as u32;
                let k = self.index.compound_pop(pivot_doc) as u32;

                // Exact Tanimoto upper bound from the known intersection size.
                let denom = p + k - c;
                let ub = if denom == 0 { 0.0 } else { c as f32 / denom as f32 };

                if ub >= threshold {
                    if let Some(fp) = get_fingerprint(pivot_doc) {
                        let (score, meets) =
                            tanimoto_with_threshold(&query.query_fp, &fp, threshold);
                        docs_evaluated += 1;
                        if meets {
                            heap.push(Candidate {
                                score: ordered_float::OrderedFloat(score),
                                doc_id: pivot_doc,
                            });
                            if heap.len() > query.top_k {
                                heap.pop();
                            }
                            if heap.len() == query.top_k {
                                if let Some(worst) = heap.peek() {
                                    threshold = threshold.max(worst.score.0);
                                }
                            }
                        }
                    }
                } else {
                    pivots_skipped += 1;
                }

                // Advance every cursor sitting on the pivot past it.
                for cursor in cursors.iter_mut() {
                    if cursor.current_doc() == Some(pivot_doc) {
                        cursor.advance_to(pivot_doc + 1);
                    }
                }
            } else {
                // Not aligned: skip a trailing cursor forward to the pivot. Docs
                // below the pivot appear in < θ lists, so none can qualify.
                pivots_skipped += 1;
                for cursor in cursors.iter_mut() {
                    if cursor.current_doc().map(|d| d < pivot_doc).unwrap_or(false) {
                        cursor.advance_to(pivot_doc);
                        break;
                    }
                }
            }
        }

        debug!(
            "WAND stats: evaluated={} pivots_skipped={} results={}",
            docs_evaluated, pivots_skipped, heap.len()
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

    /// Exhaustive Tanimoto scan: ground truth for the WAND pivot algorithm.
    ///
    /// Requires score > 0: a doc sharing zero query bits appears in none of the
    /// query's posting lists, so no inverted-index search can enumerate it. This
    /// matches WAND semantics — only docs with positive overlap are candidates.
    fn brute_force(
        query_fp: &Fingerprint,
        fps: &[Fingerprint],
        threshold: f32,
    ) -> Vec<(u32, f32)> {
        use crate::search::tanimoto::tanimoto;
        let mut out: Vec<(u32, f32)> = fps
            .iter()
            .enumerate()
            .map(|(i, fp)| (i as u32, tanimoto(query_fp, fp)))
            .filter(|(_, s)| *s >= threshold && *s > 0.0)
            .collect();
        out.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        out
    }

    #[test]
    fn test_wand_matches_brute_force() {
        // A chemically diverse set so query bits range from rare to common.
        let smiles = [
            "CCO", "CCCO", "CCCCO", "CCCCCO", "c1ccccc1", "Cc1ccccc1",
            "CN", "CC(=O)O", "CNC(=O)c1ccccc1", "c1ccncc1", "OCC(O)CO",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CC(=O)Oc1ccccc1C(=O)O",
            "C1CCCCC1", "c1ccc2ccccc2c1", "NCCO", "CCN(CC)CC",
        ];
        let (fps, index, _tmp) = build_index_from_smiles(&smiles);
        let engine = BmwEngine::new(&index);

        // Query with each compound, at several thresholds, and compare the
        // set of (doc_id, score) above threshold against the exhaustive scan.
        for q in &["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "NCCO", "CCN(CC)CC"] {
            let query_fp = compute_morgan_fp(q);
            for &t in &[0.0f32, 0.2, 0.5, 0.75, 0.9] {
                let expected = brute_force(&query_fp, &fps, t);
                let query = SimilarityQuery::new(query_fp, t, fps.len());
                let got = engine
                    .search(&query, |doc_id| fps.get(doc_id as usize).copied())
                    .unwrap();

                assert_eq!(
                    got.len(),
                    expected.len(),
                    "count mismatch for query {} @ t={}: wand={} brute={}",
                    q, t, got.len(), expected.len()
                );

                // Compare as doc_id → score maps (ordering of ties may differ).
                let mut got_sorted = got.clone();
                got_sorted.sort_by_key(|(d, _)| *d);
                let mut exp_sorted = expected.clone();
                exp_sorted.sort_by_key(|(d, _)| *d);
                for ((gd, gs), (ed, es)) in got_sorted.iter().zip(exp_sorted.iter()) {
                    assert_eq!(gd, ed, "doc_id mismatch for query {} @ t={}", q, t);
                    assert!(
                        (gs - es).abs() < 1e-6,
                        "score mismatch for doc {} query {} @ t={}: {} vs {}",
                        gd, q, t, gs, es
                    );
                }
            }
        }
    }

    #[test]
    fn test_wand_topk_truncation() {
        let smiles = [
            "CCO", "CCCO", "CCCCO", "CCCCCO", "CCCCCCO", "CCCCCCCO", "OCCO",
        ];
        let (fps, index, _tmp) = build_index_from_smiles(&smiles);
        let engine = BmwEngine::new(&index);
        let query_fp = compute_morgan_fp("CCCO");

        // top_k smaller than the number of qualifying hits must return exactly
        // the top_k highest scores from the exhaustive ranking.
        let expected = brute_force(&query_fp, &fps, 0.0);
        for &k in &[1usize, 2, 3] {
            let query = SimilarityQuery::new(query_fp, 0.0, k);
            let got = engine
                .search(&query, |doc_id| fps.get(doc_id as usize).copied())
                .unwrap();
            assert_eq!(got.len(), k);
            // The k-th best score from WAND must equal the k-th best brute score.
            let got_scores: Vec<f32> = got.iter().map(|(_, s)| *s).collect();
            for (i, gs) in got_scores.iter().enumerate() {
                assert!(
                    (gs - expected[i].1).abs() < 1e-6,
                    "rank {} score mismatch k={}: {} vs {}",
                    i, k, gs, expected[i].1
                );
            }
        }
    }
}
