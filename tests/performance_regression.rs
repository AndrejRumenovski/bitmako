//! Performance-regression guard: asserts algorithmic/pruning-ratio
//! invariants rather than wall-clock timing, which would be flaky across
//! machines and CI runners. A future change that accidentally degrades
//! WAND's pruning (e.g. a broken pivot calculation that falls back to
//! evaluating everything) should fail this test long before anyone notices
//! a slowdown on the real 1.36B-compound corpus.

mod common;

use bitmako::search::query::SimilarityQuery;
use common::{build_corpus_from_fingerprints, random_fingerprint, Xorshift64};

const CORPUS_SIZE: usize = 4000;
/// Bits set per synthetic fingerprint — chosen to roughly mirror ECFP4's
/// typical sparsity (tens of bits out of 1024), not aiming for chemical
/// realism, just a nontrivial, reproducible degree of overlap between docs.
const BITS_PER_FP: usize = 40;

fn synthetic_corpus() -> common::Corpus {
    let mut rng = Xorshift64::new(0xC0FFEE_u64);
    let fps = (0..CORPUS_SIZE).map(|_| random_fingerprint(&mut rng, BITS_PER_FP)).collect();
    build_corpus_from_fingerprints(fps)
}

#[test]
fn selective_high_threshold_query_evaluates_a_small_fraction_of_the_corpus() {
    let corpus = synthetic_corpus();
    let mut rng = Xorshift64::new(0xA5A5A5_u64);
    let query_fp = random_fingerprint(&mut rng, BITS_PER_FP);

    // At a high threshold against uniformly-random sparse fingerprints, the
    // expected overlap between any two docs is tiny (~40*40/1024 ≈ 1.6
    // shared bits), so almost nothing can reach a 0.5 Tanimoto threshold.
    // WAND should therefore evaluate only a small fraction of the corpus —
    // if a future change broke pivot pruning and fell back to a near-linear
    // scan, this bound would catch it.
    let query = SimilarityQuery::new(query_fp, 0.5, 10);
    let (_, stats) = corpus.searcher.search_with_stats(&query).unwrap();

    assert_eq!(stats.corpus_size, CORPUS_SIZE as u64);
    assert!(
        stats.eval_fraction() < 0.05,
        "WAND evaluated {:.2}% of the corpus for a selective query — pruning may be broken \
         (evaluated={}, corpus={})",
        stats.eval_fraction() * 100.0,
        stats.docs_evaluated,
        stats.corpus_size
    );
}

#[test]
fn raising_the_threshold_never_increases_result_count() {
    // With top_k large enough to never truncate, a stricter (higher)
    // Tanimoto threshold can only keep the same or fewer results than a
    // looser one — a correctness invariant independent of pruning strategy.
    let corpus = synthetic_corpus();
    let mut rng = Xorshift64::new(0xBEEF01_u64);
    let query_fp = random_fingerprint(&mut rng, BITS_PER_FP);

    let mut prev_count = usize::MAX;
    for &t in &[0.05f32, 0.15, 0.3, 0.5, 0.75] {
        let query = SimilarityQuery::new(query_fp, t, CORPUS_SIZE);
        let results = corpus.searcher.search(&query).unwrap();
        assert!(
            results.len() <= prev_count,
            "result count rose from {prev_count} to {} when threshold increased to {t}",
            results.len()
        );
        prev_count = results.len();
    }
}

#[test]
fn zero_threshold_full_recall_still_terminates_and_matches_corpus_bound() {
    // threshold=0.0 is the least selective case (WAND can't prune much) —
    // docs_evaluated must never exceed the corpus size, and the search must
    // still terminate promptly (this test itself is the timeout: if the
    // early-exit or pivot logic regressed into an infinite loop, `cargo
    // test` would hang here rather than silently pass).
    let corpus = synthetic_corpus();
    let mut rng = Xorshift64::new(0xDEAD10CC_u64);
    let query_fp = random_fingerprint(&mut rng, BITS_PER_FP);

    let query = SimilarityQuery::new(query_fp, 0.0, CORPUS_SIZE);
    let (_, stats) = corpus.searcher.search_with_stats(&query).unwrap();
    assert!(stats.docs_evaluated <= CORPUS_SIZE as u64);
}
