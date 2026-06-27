//! Similarity and property search engine.

pub mod query;
pub mod tanimoto;
pub mod wand;

use std::path::Path;

use tracing::info;

use crate::error::Result;
use crate::etl::fingerprint::{compute_morgan_fp, Fingerprint};
use crate::index::IndexReader;
use crate::search::query::SimilarityQuery;
use crate::search::wand::BmwEngine;

/// High-level search interface backed by the BMW engine and a Lance dataset.
///
/// For large-scale search, fingerprints are fetched lazily from the Lance dataset
/// during BMW evaluation (only accessed when the engine decides a document is a
/// viable candidate above the current threshold).
pub struct Searcher {
    index: IndexReader,
    /// Compact flat fingerprint store for sub-linear access by doc_id
    fp_store: Vec<Fingerprint>,
}

impl Searcher {
    /// Load the index from disk and use the provided in-memory fingerprint store.
    pub fn open(index_path: &Path, fp_store: Vec<Fingerprint>) -> Result<Self> {
        let index = IndexReader::open(index_path)?;
        info!("Searcher ready: {} compounds indexed", index.num_compounds);
        Ok(Searcher { index, fp_store })
    }

    /// Construct a Searcher from an already-loaded IndexReader.
    pub fn open_from_index(index: IndexReader, fp_store: Vec<Fingerprint>) -> Self {
        info!("Searcher ready: {} compounds indexed", index.num_compounds);
        Searcher { index, fp_store }
    }

    /// Execute a similarity search returning top-k results.
    pub fn search(&self, query: &SimilarityQuery) -> Result<Vec<(u32, f32)>> {
        let engine = BmwEngine::new(&self.index);
        let results = engine.search(query, |doc_id| {
            self.fp_store.get(doc_id as usize).copied()
        })?;
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
