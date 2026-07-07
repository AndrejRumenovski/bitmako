//! Similarity analysis: turns a raw Tanimoto score into an explainable
//! bit-level breakdown of *why* two fingerprints matched, for display
//! alongside search results (the "Similarity Analysis" panel).
//!
//! This is a pure post-processing step over two already-fetched fingerprints —
//! it does not touch the WAND pivot loop, the posting lists, or block-max
//! pruning, so it has no effect on search performance. It only ever runs over
//! the (small) top-k result set a search already returned, never the corpus:
//! computing it for 50 results costs 50 extra O(1) mmap reads plus 50 cheap
//! 16-word bit scans, negligible next to the search itself.

use crate::etl::fingerprint::{Fingerprint, FP_WORDS};

/// A bit-level breakdown of how a candidate fingerprint compares to a query
/// fingerprint, plus a human-readable explanation of the match.
#[derive(Debug, Clone, PartialEq)]
pub struct SimilarityAnalysis {
    /// Exact Tanimoto similarity (`shared_bits / total_bits`). Numerically
    /// identical to the score WAND already computed for this candidate —
    /// same fingerprints, same math.
    pub tanimoto: f32,
    /// Bits set in both the query and the candidate fingerprint.
    pub shared_bits: u32,
    /// Bits set in the query fingerprint but not in the candidate's.
    pub query_unique_bits: u32,
    /// Bits set in the candidate fingerprint but not in the query's.
    pub candidate_unique_bits: u32,
    /// 1–2 sentence natural-language description of the match, generated
    /// directly from the bit counts above — not a canned per-molecule string.
    pub explanation: String,
}

impl SimilarityAnalysis {
    /// Total number of distinct bits set across either fingerprint
    /// (`shared + query_unique + candidate_unique`) — the Tanimoto union.
    pub fn total_bits(&self) -> u32 {
        self.shared_bits + self.query_unique_bits + self.candidate_unique_bits
    }
}

/// Compare a query fingerprint against a candidate's, computing the exact
/// bit-level breakdown and a generated explanation.
///
/// Cost is O(`FP_WORDS`) = 16 word ops — the same inner loop
/// [`tanimoto`](crate::search::tanimoto::tanimoto) already uses.
pub fn analyze(query: &Fingerprint, candidate: &Fingerprint) -> SimilarityAnalysis {
    let mut shared = 0u32;
    let mut query_only = 0u32;
    let mut candidate_only = 0u32;

    for i in 0..FP_WORDS {
        shared += (query[i] & candidate[i]).count_ones();
        query_only += (query[i] & !candidate[i]).count_ones();
        candidate_only += (candidate[i] & !query[i]).count_ones();
    }

    let union = shared + query_only + candidate_only;
    let tanimoto = if union == 0 { 0.0 } else { shared as f32 / union as f32 };

    SimilarityAnalysis {
        tanimoto,
        shared_bits: shared,
        query_unique_bits: query_only,
        candidate_unique_bits: candidate_only,
        explanation: explain(tanimoto, shared, query_only, candidate_only),
    }
}

/// Qualitative word for the explanation sentence. This banding is a display
/// heuristic only — the numbers behind it (in the sentence and the panel's
/// other fields) are always the real, freshly computed values, not looked up
/// from a table.
fn similarity_tier(tanimoto: f32) -> &'static str {
    match tanimoto {
        t if t >= 0.7 => "high",
        t if t >= 0.4 => "moderate",
        t if t >= 0.2 => "low but notable",
        _ => "minimal",
    }
}

/// Build the 1–2 sentence explanation from the real bit counts for this pair.
fn explain(tanimoto: f32, shared: u32, query_only: u32, candidate_only: u32) -> String {
    let pct = (tanimoto * 100.0).round() as i64;
    let tier = similarity_tier(tanimoto);
    let total = shared + query_only + candidate_only;

    format!(
        "This molecule shares {pct}% of its structural fingerprint with the query \
         ({shared} of {total} total bits), indicating a {tier} degree of structural \
         similarity. {query_only} bit{qs} appear only in the query and {candidate_only} \
         bit{cs} appear only in this molecule.",
        qs = if query_only == 1 { "" } else { "s" },
        cs = if candidate_only == 1 { "" } else { "s" },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::fingerprint::compute_morgan_fp;
    use crate::search::tanimoto::tanimoto;

    #[test]
    fn identical_fingerprints_are_fully_shared() {
        let fp = compute_morgan_fp("CC(=O)Oc1ccccc1C(=O)O");
        let a = analyze(&fp, &fp);
        assert_eq!(a.query_unique_bits, 0);
        assert_eq!(a.candidate_unique_bits, 0);
        assert!((a.tanimoto - 1.0).abs() < 1e-6);
        assert!(a.explanation.contains("100%"));
    }

    #[test]
    fn zero_fingerprints_do_not_panic() {
        let zero = [0u64; 16];
        let a = analyze(&zero, &zero);
        assert_eq!(a.tanimoto, 0.0);
        assert_eq!(a.shared_bits, 0);
        assert_eq!(a.total_bits(), 0);
        assert!(a.explanation.contains("0%"));
    }

    #[test]
    fn disjoint_fingerprints_have_no_shared_bits() {
        let mut a = [0u64; 16];
        let mut b = [0u64; 16];
        a[0] = 0x00FF_00FF_00FF_00FFu64;
        b[0] = 0xFF00_FF00_FF00_FF00u64;
        let analysis = analyze(&a, &b);
        assert_eq!(analysis.shared_bits, 0);
        assert_eq!(analysis.query_unique_bits, 32);
        assert_eq!(analysis.candidate_unique_bits, 32);
        assert_eq!(analysis.tanimoto, 0.0);
    }

    #[test]
    fn bit_counts_match_manual_partial_overlap() {
        let mut a = [0u64; 16];
        let mut b = [0u64; 16];
        a[0] = 0b0000_1111u64; // query bits 0-3
        b[0] = 0b0011_1100u64; // candidate bits 2-5
        let analysis = analyze(&a, &b);
        assert_eq!(analysis.shared_bits, 2); // bits 2,3
        assert_eq!(analysis.query_unique_bits, 2); // bits 0,1
        assert_eq!(analysis.candidate_unique_bits, 2); // bits 4,5
        assert_eq!(analysis.total_bits(), 6);
        assert!((analysis.tanimoto - 2.0 / 6.0).abs() < 1e-6);
    }

    #[test]
    fn tanimoto_matches_the_shared_tanimoto_function() {
        // The analysis module must never disagree with the score WAND already
        // returned for the same pair — it's the same comparison, computed twice.
        let pairs = [
            ("CCO", "CCCO"),
            ("c1ccccc1", "Cc1ccccc1"),
            ("CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1"),
        ];
        for (q, c) in pairs {
            let qfp = compute_morgan_fp(q);
            let cfp = compute_morgan_fp(c);
            let a = analyze(&qfp, &cfp);
            assert!((a.tanimoto - tanimoto(&qfp, &cfp)).abs() < 1e-6);
        }
    }

    #[test]
    fn explanation_reflects_similarity_tier() {
        let fp = compute_morgan_fp("CC(=O)Oc1ccccc1C(=O)O");
        let high = analyze(&fp, &fp);
        assert!(high.explanation.contains("high degree"));

        let mut disjoint_b = [0u64; 16];
        disjoint_b[15] = 1u64 << 63; // a bit almost certainly unset in `fp`
        let minimal = analyze(&fp, &disjoint_b);
        assert!(minimal.explanation.contains("minimal degree") || minimal.tanimoto < 0.2);
    }
}
