//! Tanimoto (Jaccard) similarity between 1024-bit binary fingerprints.
//!
//! Tanimoto(A, B) = |A ∩ B| / |A ∪ B|
//!                = popcount(A AND B) / popcount(A OR B)
//!
//! LLVM will autovectorize the inner loops over the 16 u64 words
//! when compiled with RUSTFLAGS="-C target-cpu=native".

use crate::etl::fingerprint::{Fingerprint, FP_WORDS};

/// Compute Tanimoto similarity between two 1024-bit fingerprints.
/// Returns a value in [0.0, 1.0].
#[inline]
pub fn tanimoto(a: &Fingerprint, b: &Fingerprint) -> f32 {
    let mut intersection = 0u32;
    let mut union_bits = 0u32;

    // LLVM autovectorizes this loop: each .count_ones() maps to POPCNT instruction
    for i in 0..FP_WORDS {
        intersection += (a[i] & b[i]).count_ones();
        union_bits += (a[i] | b[i]).count_ones();
    }

    if union_bits == 0 {
        return 0.0;
    }

    intersection as f32 / union_bits as f32
}

/// Compute Tanimoto from raw word slices (avoids double-dereference overhead in hot paths)
#[inline]
pub fn tanimoto_slice(a: &[u64; FP_WORDS], b: &[u64; FP_WORDS]) -> f32 {
    tanimoto(a, b)
}

/// Check if Tanimoto meets threshold without computing exact score (early termination).
/// Returns (similarity, met_threshold).
#[inline]
pub fn tanimoto_with_threshold(a: &Fingerprint, b: &Fingerprint, threshold: f32) -> (f32, bool) {
    let mut intersection = 0u32;
    let mut union_bits = 0u32;

    // Compute exact Tanimoto — POPCNT is fast enough that branching per-word is slower
    for i in 0..FP_WORDS {
        intersection += (a[i] & b[i]).count_ones();
        union_bits += (a[i] | b[i]).count_ones();
    }

    if union_bits == 0 {
        return (0.0, false);
    }

    let sim = intersection as f32 / union_bits as f32;
    (sim, sim >= threshold)
}

/// Upper bound on Tanimoto when we know the query and candidate popcounts.
/// Useful for WAND-style pruning before computing exact Tanimoto.
///
/// Tanimoto UB = min(pop_a, pop_b) / max(pop_a, pop_b)
/// (achieved when one fingerprint is a subset of the other)
#[inline]
pub fn tanimoto_upper_bound(pop_a: u32, pop_b: u32) -> f32 {
    if pop_a == 0 || pop_b == 0 {
        return 0.0;
    }
    let min_pop = pop_a.min(pop_b);
    let max_pop = pop_a.max(pop_b);
    min_pop as f32 / max_pop as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_fingerprints() {
        let fp: Fingerprint = [0xDEADBEEFCAFEBABEu64; 16];
        assert!((tanimoto(&fp, &fp) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_zero_fingerprints() {
        let fp: Fingerprint = [0u64; 16];
        assert_eq!(tanimoto(&fp, &fp), 0.0);
    }

    #[test]
    fn test_disjoint_fingerprints() {
        let mut a: Fingerprint = [0u64; 16];
        let mut b: Fingerprint = [0u64; 16];
        a[0] = 0x00FF_00FF_00FF_00FFu64;
        b[0] = 0xFF00_FF00_FF00_FF00u64;
        // intersection = 0, union = all bits in a and b
        assert!((tanimoto(&a, &b) - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_partial_overlap() {
        let mut a: Fingerprint = [0u64; 16];
        let mut b: Fingerprint = [0u64; 16];
        // a has bits 0-3 set, b has bits 2-5 set in word 0
        a[0] = 0b0000_1111u64; // 4 bits
        b[0] = 0b0011_1100u64; // 4 bits
        // intersection = bits 2-3 = 2 bits, union = bits 0-5 = 6 bits
        let sim = tanimoto(&a, &b);
        assert!((sim - 2.0 / 6.0).abs() < 1e-6);
    }

    #[test]
    fn test_upper_bound() {
        assert!((tanimoto_upper_bound(10, 100) - 0.1).abs() < 1e-6);
        assert!((tanimoto_upper_bound(50, 50) - 1.0).abs() < 1e-6);
        assert_eq!(tanimoto_upper_bound(0, 50), 0.0);
    }
}
