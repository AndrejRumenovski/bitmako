//! Similarity and property search engine.

pub mod fp_store;
pub mod prop_store;
pub mod query;
pub mod tanimoto;
pub mod wand;

use std::path::Path;

use tracing::info;

use crate::error::Result;
use crate::etl::fingerprint::compute_morgan_fp;
use crate::index::skip::SkipIndex;
use crate::index::IndexReader;
use crate::search::fp_store::FpStore;
use crate::search::prop_store::PropStore;
use crate::search::query::SimilarityQuery;
use crate::search::wand::{BmwEngine, SearchStats};

/// High-level search interface backed by the BMW engine and a flat fingerprint store.
///
/// Fingerprints are fetched lazily from a memory-mapped flat file during BMW
/// evaluation — only accessed when the engine decides a document is a viable
/// candidate above the current threshold, so the resident set stays small.
pub struct Searcher {
    index: IndexReader,
    /// Skip index for streaming, block-at-a-time posting-list traversal.
    skip: SkipIndex,
    /// Memory-mapped flat fingerprint store for O(1) access by doc_id
    fp_store: FpStore,
    /// Optional memory-mapped flat property store. When present and the query
    /// carries property filters, BMW screens properties inside the pivot loop —
    /// before paying for the fingerprint fetch — so no over-fetch is needed.
    prop_store: Option<PropStore>,
}

impl Searcher {
    /// Load the index, skip index, and fingerprint store from disk.
    pub fn open(index_path: &Path, skip_path: &Path, fp_store_path: &Path) -> Result<Self> {
        let index = IndexReader::open(index_path)?;
        let skip = SkipIndex::open(skip_path)?;
        let fp_store = FpStore::open(fp_store_path)?;
        info!(
            "Searcher ready: {} compounds indexed, {} fingerprints in store",
            index.num_compounds,
            fp_store.len()
        );
        Ok(Searcher { index, skip, fp_store, prop_store: None })
    }

    /// Construct a Searcher from already-loaded components.
    pub fn open_from_index(index: IndexReader, skip: SkipIndex, fp_store: FpStore) -> Self {
        info!(
            "Searcher ready: {} compounds indexed, {} fingerprints in store",
            index.num_compounds,
            fp_store.len()
        );
        Searcher { index, skip, fp_store, prop_store: None }
    }

    /// Attach a flat property store, enabling in-loop property pre-filtering.
    pub fn with_prop_store(mut self, prop_store: PropStore) -> Self {
        info!("Property store attached: {} records", prop_store.len());
        self.prop_store = Some(prop_store);
        self
    }

    /// Number of compounds indexed.
    pub fn num_compounds(&self) -> u32 {
        self.index.num_compounds
    }

    /// Whether a property store is attached, enabling `--mw-max`/`--logp-max`
    /// filters to be screened in-loop instead of silently ignored.
    pub fn has_prop_store(&self) -> bool {
        self.prop_store.is_some()
    }

    /// Execute a similarity search returning top-k results.
    ///
    /// When a property store is attached and the query carries property filters,
    /// the filter is evaluated inside the BMW pivot loop as a cheap pre-screen, so
    /// only compounds satisfying both the Tanimoto threshold *and* the property
    /// filters reach the heap — `top_k` results come back already filtered, with
    /// no over-fetch.
    pub fn search(&self, query: &SimilarityQuery) -> Result<Vec<(u32, f32)>> {
        let engine = BmwEngine::new(&self.index, &self.skip);
        let results = match &self.prop_store {
            Some(prop_store) if !query.property_filters.is_empty() => engine.search_filtered(
                query,
                |doc_id| self.fp_store.get(doc_id),
                |doc_id| {
                    prop_store
                        .get(doc_id)
                        .map(|props| query.filter_passes(&props))
                        .unwrap_or(false)
                },
            )?,
            _ => engine.search(query, |doc_id| self.fp_store.get(doc_id))?,
        };
        Ok(results)
    }

    /// Like [`search`] but also returns pruning diagnostics.
    pub fn search_with_stats(&self, query: &SimilarityQuery) -> Result<(Vec<(u32, f32)>, SearchStats)> {
        let engine = BmwEngine::new(&self.index, &self.skip);
        match &self.prop_store {
            Some(prop_store) if !query.property_filters.is_empty() => engine.search_filtered_with_stats(
                query,
                |doc_id| self.fp_store.get(doc_id),
                |doc_id| {
                    prop_store
                        .get(doc_id)
                        .map(|props| query.filter_passes(&props))
                        .unwrap_or(false)
                },
            ),
            _ => engine.search_with_stats(query, |doc_id| self.fp_store.get(doc_id)),
        }
    }

    /// Convenience: search by SMILES string instead of pre-computed fingerprint.
    pub fn search_by_smiles(
        &self,
        smiles: &str,
        threshold: f32,
        top_k: usize,
    ) -> Result<Vec<(u32, f32)>> {
        let query_fp = compute_morgan_fp(smiles);
        let query = SimilarityQuery::new(query_fp, threshold, top_k);
        self.search(&query)
    }
}
