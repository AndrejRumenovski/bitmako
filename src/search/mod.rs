//! Similarity and property search engine.

pub mod fp_store;
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
use crate::search::query::SimilarityQuery;
use crate::search::wand::BmwEngine;

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
        Ok(Searcher { index, skip, fp_store })
    }

    /// Construct a Searcher from already-loaded components.
    pub fn open_from_index(index: IndexReader, skip: SkipIndex, fp_store: FpStore) -> Self {
        info!(
            "Searcher ready: {} compounds indexed, {} fingerprints in store",
            index.num_compounds,
            fp_store.len()
        );
        Searcher { index, skip, fp_store }
    }

    /// Execute a similarity search returning top-k results.
    pub fn search(&self, query: &SimilarityQuery) -> Result<Vec<(u32, f32)>> {
        let engine = BmwEngine::new(&self.index, &self.skip);
        let results = engine.search(query, |doc_id| self.fp_store.get(doc_id))?;
        Ok(results)
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
