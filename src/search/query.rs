//! Query types for the BMW similarity search engine.
//!
//! Supports:
//!   - Structural similarity queries (Tanimoto threshold on Morgan fingerprints)
//!   - Scalar property filters (LogP, MW, RotBonds, HeavyAtoms)
//!   - Conjunctive combinations of the above

use crate::etl::fingerprint::{fp_popcount, Fingerprint};
use crate::etl::properties::MolecularProperties;
use crate::error::{BitMakoError, Result};

/// A range filter on a scalar molecular property.
#[derive(Debug, Clone)]
pub struct PropertyFilter {
    pub field: PropertyField,
    pub min: Option<f32>,
    pub max: Option<f32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropertyField {
    MolWeight,
    LogP,
    RotatableBonds,
    HeavyAtomCount,
    RingCount,
}

impl PropertyFilter {
    pub fn passes(&self, props: &MolecularProperties) -> bool {
        let value = match self.field {
            PropertyField::MolWeight => props.mw,
            PropertyField::LogP => props.logp,
            PropertyField::RotatableBonds => props.rot_bonds as f32,
            PropertyField::HeavyAtomCount => props.heavy_atoms as f32,
            PropertyField::RingCount => props.ring_count as f32,
        };
        self.min.map(|m| value >= m).unwrap_or(true)
            && self.max.map(|m| value <= m).unwrap_or(true)
    }
}

/// A complete similarity + property search query
#[derive(Debug, Clone)]
pub struct SimilarityQuery {
    /// Query fingerprint (1024-bit Morgan FP of the probe compound)
    pub query_fp: Fingerprint,
    /// Pre-computed popcount of the query fingerprint
    pub query_pop: u32,
    /// Minimum Tanimoto similarity threshold (0.0..=1.0)
    pub threshold: f32,
    /// Maximum number of results to return (top-k)
    pub top_k: usize,
    /// Optional conjunctive property filters applied after structural screening
    pub property_filters: Vec<PropertyFilter>,
}

impl SimilarityQuery {
    pub fn new(query_fp: Fingerprint, threshold: f32, top_k: usize) -> Self {
        let query_pop = fp_popcount(&query_fp);
        SimilarityQuery {
            query_fp,
            query_pop,
            threshold,
            top_k,
            property_filters: Vec::new(),
        }
    }

    /// Builder: add a property filter
    pub fn with_filter(mut self, filter: PropertyFilter) -> Self {
        self.property_filters.push(filter);
        self
    }

    /// Builder: add `mw_max`/`logp_max` upper-bound filters when present. This is
    /// the common case for CLI/API callers exposing `--mw-max`/`--logp-max` (or
    /// their JSON equivalents) as optional, independent filters.
    pub fn with_mw_logp_max(mut self, mw_max: Option<f32>, logp_max: Option<f32>) -> Self {
        if let Some(max) = mw_max {
            self = self.with_filter(PropertyFilter { field: PropertyField::MolWeight, min: None, max: Some(max) });
        }
        if let Some(max) = logp_max {
            self = self.with_filter(PropertyFilter { field: PropertyField::LogP, min: None, max: Some(max) });
        }
        self
    }

    /// Builder: add a Lipinski-style drug-like filter
    pub fn with_lipinski_filter(self) -> Self {
        self
            .with_filter(PropertyFilter { field: PropertyField::MolWeight, min: None, max: Some(500.0) })
            .with_filter(PropertyFilter { field: PropertyField::LogP, min: None, max: Some(5.0) })
            .with_filter(PropertyFilter { field: PropertyField::RotatableBonds, min: None, max: Some(10.0) })
    }

    /// Returns the set of query bit indices (bits set to 1 in the query FP).
    /// These are the active posting lists to scan.
    pub fn active_bits(&self) -> Vec<usize> {
        let mut bits = Vec::new();
        for (word_idx, &word) in self.query_fp.iter().enumerate() {
            let mut w = word;
            while w != 0 {
                let trailing = w.trailing_zeros() as usize;
                bits.push(word_idx * 64 + trailing);
                w &= w - 1;
            }
        }
        bits
    }

    /// True if properties satisfy all filters
    pub fn filter_passes(&self, props: &MolecularProperties) -> bool {
        self.property_filters.iter().all(|f| f.passes(props))
    }

    /// Tanimoto upper bound given a candidate popcount.
    /// Used for WAND-style early termination.
    #[inline]
    pub fn tanimoto_upper_bound(&self, candidate_pop: u32) -> f32 {
        crate::search::tanimoto::tanimoto_upper_bound(self.query_pop, candidate_pop)
    }
}

/// Validate query parameters
pub fn validate_query(query: &SimilarityQuery) -> Result<()> {
    if !(0.0..=1.0).contains(&query.threshold) {
        return Err(BitMakoError::Query(format!(
            "Tanimoto threshold {} is out of range [0.0, 1.0]",
            query.threshold
        )));
    }
    if query.top_k == 0 {
        return Err(BitMakoError::Query("top_k must be >= 1".into()));
    }
    if query.query_pop == 0 {
        return Err(BitMakoError::Query(
            "Query fingerprint is all zeros — no results possible".into()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::fingerprint::compute_morgan_fp;

    #[test]
    fn test_active_bits() {
        let mut fp = [0u64; 16];
        fp[0] = 0b1011u64; // bits 0, 1, 3
        let query = SimilarityQuery::new(fp, 0.7, 10);
        let bits = query.active_bits();
        assert!(bits.contains(&0));
        assert!(bits.contains(&1));
        assert!(bits.contains(&3));
        assert!(!bits.contains(&2));
    }

    #[test]
    fn test_lipinski_filter() {
        use crate::etl::properties::MolecularProperties;
        let fp = compute_morgan_fp("CCO");
        let query = SimilarityQuery::new(fp, 0.7, 10).with_lipinski_filter();
        let props = MolecularProperties { mw: 300.0, logp: 2.0, rot_bonds: 3, heavy_atoms: 20, ring_count: 1 };
        assert!(query.filter_passes(&props));

        let heavy = MolecularProperties { mw: 600.0, ..props };
        assert!(!query.filter_passes(&heavy));
    }

    #[test]
    fn test_validate_bad_threshold() {
        let fp = compute_morgan_fp("CCO");
        let query = SimilarityQuery::new(fp, 1.5, 10);
        assert!(validate_query(&query).is_err());
    }
}
