//! Shared Lance dataset lookup: resolving doc_ids to SMILES/properties.
//!
//! Every caller that needs to show a human a search result — single-query
//! `search --lance`, batch `search-batch --lance`, and the HTTP API's `/search` —
//! does the same thing: take the WAND-ranked doc_ids as row indices into the Lance
//! dataset and pull out `compound_id`/`smiles`/`mw`/`logp`/`rot_bonds`/`heavy_atoms`/
//! `ring_count`. This was previously reimplemented at each call site; this module
//! is the one place that does it.

use crate::error::{LanceResultExt, Result};
use crate::etl::properties::MolecularProperties;

/// One resolved row: a compound's catalog ID, SMILES, and molecular properties.
#[derive(Debug, Clone)]
pub struct ResolvedCompound {
    pub compound_id: String,
    pub smiles: String,
    pub properties: MolecularProperties,
}

const PROPERTY_COLUMNS: &[&str] =
    &["compound_id", "smiles", "mw", "logp", "rot_bonds", "heavy_atoms", "ring_count"];

/// Fetch compound_id/SMILES/properties for each doc_id in `doc_ids`, in the same
/// order, via a single `Dataset::take` call (random-access row lookup, no scan).
pub async fn resolve_compounds(
    dataset: &lance::dataset::Dataset,
    doc_ids: &[u32],
) -> Result<Vec<ResolvedCompound>> {
    use arrow_array::cast::AsArray;
    use arrow_array::types::{Float32Type, UInt32Type};

    if doc_ids.is_empty() {
        return Ok(Vec::new());
    }

    let row_indices: Vec<u64> = doc_ids.iter().map(|&d| d as u64).collect();
    let projection = dataset.schema().project(PROPERTY_COLUMNS).lance_err()?;
    let batch = dataset.take(&row_indices, projection).await.lance_err()?;

    let cid_col = batch.column_by_name("compound_id").unwrap().as_string::<i32>();
    let smi_col = batch.column_by_name("smiles").unwrap().as_string::<i32>();
    let mw_col = batch.column_by_name("mw").unwrap().as_primitive::<Float32Type>();
    let logp_col = batch.column_by_name("logp").unwrap().as_primitive::<Float32Type>();
    let rot_col = batch.column_by_name("rot_bonds").unwrap().as_primitive::<UInt32Type>();
    let heavy_col = batch.column_by_name("heavy_atoms").unwrap().as_primitive::<UInt32Type>();
    let ring_col = batch.column_by_name("ring_count").unwrap().as_primitive::<UInt32Type>();

    Ok((0..doc_ids.len())
        .map(|i| ResolvedCompound {
            compound_id: cid_col.value(i).to_string(),
            smiles: smi_col.value(i).to_string(),
            properties: MolecularProperties {
                mw: mw_col.value(i),
                logp: logp_col.value(i),
                rot_bonds: rot_col.value(i),
                heavy_atoms: heavy_col.value(i),
                ring_count: ring_col.value(i),
            },
        })
        .collect())
}
