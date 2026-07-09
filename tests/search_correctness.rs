//! Black-box search-correctness tests over the public `Searcher` API.
//!
//! These cross-check the full WAND pipeline (index + skip index + fp store +
//! prop store, wired together exactly as `Searcher::open` does for the real
//! 1.4B-compound corpus) against an exhaustive brute-force Tanimoto scan, on
//! synthetic corpora built purely in-memory/temp-files — no dependency on the
//! real `data/` files.

mod common;

use bitmako::etl::fingerprint::compute_morgan_fp;
use bitmako::search::query::{PropertyField, PropertyFilter, SimilarityQuery};
use common::{brute_force, build_corpus, DIVERSE_SMILES};

#[test]
fn wand_matches_brute_force_across_queries_and_thresholds() {
    let corpus = build_corpus(DIVERSE_SMILES);

    for q in DIVERSE_SMILES {
        let query_fp = compute_morgan_fp(q);
        for &t in &[0.0f32, 0.1, 0.2, 0.5, 0.75, 0.9, 1.0] {
            let expected = brute_force(&query_fp, &corpus.fps, t);
            let query = SimilarityQuery::new(query_fp, t, corpus.fps.len());
            let got = corpus.searcher.search(&query).expect("search must succeed");

            assert_eq!(
                got.len(),
                expected.len(),
                "count mismatch for query {q} @ t={t}: wand={} brute={}",
                got.len(),
                expected.len()
            );

            let mut got_sorted = got.clone();
            got_sorted.sort_by_key(|(d, _)| *d);
            let mut exp_sorted = expected.clone();
            exp_sorted.sort_by_key(|(d, _)| *d);
            for ((gd, gs), (ed, es)) in got_sorted.iter().zip(exp_sorted.iter()) {
                assert_eq!(gd, ed, "doc_id mismatch for query {q} @ t={t}");
                assert!(
                    (gs - es).abs() < 1e-6,
                    "score mismatch for doc {gd} query {q} @ t={t}: {gs} vs {es}"
                );
            }
        }
    }
}

#[test]
fn results_are_sorted_descending_by_score() {
    let corpus = build_corpus(DIVERSE_SMILES);
    let query_fp = compute_morgan_fp("CC(=O)Oc1ccccc1C(=O)O");
    let query = SimilarityQuery::new(query_fp, 0.0, corpus.fps.len());
    let results = corpus.searcher.search(&query).unwrap();

    for window in results.windows(2) {
        assert!(window[0].1 >= window[1].1, "results not sorted descending");
    }
}

#[test]
fn top_k_truncates_to_the_highest_scoring_subset() {
    let corpus = build_corpus(DIVERSE_SMILES);
    let query_fp = compute_morgan_fp("CCCO");
    let expected = brute_force(&query_fp, &corpus.fps, 0.0);

    for &k in &[1usize, 3, 5, 10] {
        let query = SimilarityQuery::new(query_fp, 0.0, k);
        let got = corpus.searcher.search(&query).unwrap();
        assert_eq!(got.len(), k.min(expected.len()));
        for (i, (_, gs)) in got.iter().enumerate() {
            assert!(
                (gs - expected[i].1).abs() < 1e-6,
                "rank {i} score mismatch k={k}: {gs} vs {}",
                expected[i].1
            );
        }
    }
}

#[test]
fn search_is_deterministic_across_repeated_runs() {
    let corpus = build_corpus(DIVERSE_SMILES);
    let query_fp = compute_morgan_fp("c1ccccc1");
    let query = SimilarityQuery::new(query_fp, 0.1, 10);

    let first = corpus.searcher.search(&query).unwrap();
    for _ in 0..5 {
        let again = corpus.searcher.search(&query).unwrap();
        assert_eq!(first, again, "identical query must return identical results every time");
    }
}

#[test]
fn property_filtered_search_matches_brute_force_conjunction() {
    let corpus = build_corpus(DIVERSE_SMILES);
    let query_fp = compute_morgan_fp("CCO");
    let mw_max = 100.0f32;

    let query = SimilarityQuery::new(query_fp, 0.0, corpus.fps.len()).with_filter(PropertyFilter {
        field: PropertyField::MolWeight,
        min: None,
        max: Some(mw_max),
    });
    let got = corpus.searcher.search(&query).unwrap();

    // Brute-force conjunction: Tanimoto > 0 AND mw <= mw_max.
    let mut expected: Vec<(u32, f32)> = corpus
        .fps
        .iter()
        .zip(corpus.props.iter())
        .enumerate()
        .filter_map(|(i, (fp, props))| {
            let score = bitmako::search::tanimoto::tanimoto(&query_fp, fp);
            if score > 0.0 && props.mw <= mw_max {
                Some((i as u32, score))
            } else {
                None
            }
        })
        .collect();
    expected.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    let mut got_sorted = got.clone();
    got_sorted.sort_by_key(|(d, _)| *d);
    let mut exp_sorted = expected.clone();
    exp_sorted.sort_by_key(|(d, _)| *d);
    assert_eq!(got_sorted.len(), exp_sorted.len(), "property-filtered result count mismatch");
    for ((gd, gs), (ed, es)) in got_sorted.iter().zip(exp_sorted.iter()) {
        assert_eq!(gd, ed);
        assert!((gs - es).abs() < 1e-6);
    }
}

#[test]
fn property_filter_never_returns_a_compound_that_fails_the_filter() {
    let corpus = build_corpus(DIVERSE_SMILES);
    let query_fp = compute_morgan_fp("CC(=O)Oc1ccccc1C(=O)O");
    let logp_max = 1.5f32;

    let query = SimilarityQuery::new(query_fp, 0.0, corpus.fps.len()).with_filter(PropertyFilter {
        field: PropertyField::LogP,
        min: None,
        max: Some(logp_max),
    });
    let got = corpus.searcher.search(&query).unwrap();

    for (doc_id, _) in &got {
        let props = corpus.props[*doc_id as usize];
        assert!(
            props.logp <= logp_max,
            "returned doc {doc_id} has logp={} exceeding filter max {logp_max}",
            props.logp
        );
    }
}

#[test]
fn analyze_results_bit_counts_are_consistent_with_search_scores() {
    let corpus = build_corpus(DIVERSE_SMILES);
    let query_fp = compute_morgan_fp("NCCO");
    let query = SimilarityQuery::new(query_fp, 0.0, corpus.fps.len());
    let results = corpus.searcher.search(&query).unwrap();
    let analyses = corpus.searcher.analyze_results(&query, &results);

    assert_eq!(analyses.len(), results.len(), "one analysis per result expected");
    for ((_, score), analysis) in results.iter().zip(analyses.iter()) {
        assert!(
            (analysis.tanimoto - score).abs() < 1e-6,
            "analysis tanimoto {} disagrees with search score {score}",
            analysis.tanimoto
        );
        let union = analysis.shared_bits + analysis.query_unique_bits + analysis.candidate_unique_bits;
        assert_eq!(union, analysis.total_bits());
        if union > 0 {
            let recomputed = analysis.shared_bits as f32 / union as f32;
            assert!((recomputed - analysis.tanimoto).abs() < 1e-6);
        }
    }
}

#[test]
fn exact_self_match_scores_one() {
    let corpus = build_corpus(DIVERSE_SMILES);
    let query_fp = compute_morgan_fp(DIVERSE_SMILES[0]);
    let query = SimilarityQuery::new(query_fp, 0.5, 5);
    let results = corpus.searcher.search(&query).unwrap();

    assert_eq!(results[0].0, 0, "querying a corpus member with itself should rank it first");
    assert!((results[0].1 - 1.0).abs() < 1e-5);
}

#[test]
fn empty_corpus_never_reached_but_single_compound_corpus_works() {
    let corpus = build_corpus(&["CCO"]);
    let query_fp = compute_morgan_fp("CCO");
    let query = SimilarityQuery::new(query_fp, 0.0, 10);
    let results = corpus.searcher.search(&query).unwrap();
    assert_eq!(results, vec![(0, 1.0)]);
}
