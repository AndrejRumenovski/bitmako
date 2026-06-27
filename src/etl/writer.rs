//! Arrow RecordBatch builder and Lance dataset writer.
//!
//! Converts parsed compound batches into columnar Arrow format, optionally
//! runs Vortex compression on the output arrays, then appends to a Lance dataset.
//!
//! Schema:
//!   compound_id   Utf8
//!   smiles        Utf8
//!   mw            Float32
//!   logp          Float32
//!   rot_bonds     UInt32
//!   heavy_atoms   UInt32
//!   ring_count    UInt32
//!   fingerprint   FixedSizeList(UInt64, 16)   -- 1024-bit Morgan FP

use std::path::Path;
use std::sync::Arc;

use arrow_array::{
    FixedSizeListArray, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array,
    Array,
};
use arrow_schema::{DataType, Field, Schema};
use tracing::{debug, info};

use crate::error::{BitMakoError, Result};
use crate::etl::fingerprint::FP_WORDS;
use crate::etl::parser::ParsedCompound;

/// Returns the Arrow schema for the compound dataset.
pub fn compound_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("compound_id", DataType::Utf8, false),
        Field::new("smiles", DataType::Utf8, false),
        Field::new("mw", DataType::Float32, false),
        Field::new("logp", DataType::Float32, false),
        Field::new("rot_bonds", DataType::UInt32, false),
        Field::new("heavy_atoms", DataType::UInt32, false),
        Field::new("ring_count", DataType::UInt32, false),
        Field::new(
            "fingerprint",
            DataType::FixedSizeList(
                Arc::new(Field::new("item", DataType::UInt64, false)),
                FP_WORDS as i32,
            ),
            false,
        ),
    ]))
}

/// Build an Arrow RecordBatch from a slice of ParsedCompound.
pub fn build_record_batch(compounds: &[ParsedCompound]) -> Result<RecordBatch> {
    let n = compounds.len();
    let schema = compound_schema();

    let compound_ids = StringArray::from(compounds.iter().map(|c| c.compound_id.as_str()).collect::<Vec<_>>());
    let smiles = StringArray::from(compounds.iter().map(|c| c.smiles.as_str()).collect::<Vec<_>>());
    let mw: Float32Array = compounds.iter().map(|c| c.properties.mw).collect();
    let logp: Float32Array = compounds.iter().map(|c| c.properties.logp).collect();
    let rot_bonds: UInt32Array = compounds.iter().map(|c| c.properties.rot_bonds).collect();
    let heavy_atoms: UInt32Array = compounds.iter().map(|c| c.properties.heavy_atoms).collect();
    let ring_count: UInt32Array = compounds.iter().map(|c| c.properties.ring_count).collect();

    // Flatten fingerprints into a single u64 array, then wrap as FixedSizeList
    let mut fp_values: Vec<u64> = Vec::with_capacity(n * FP_WORDS);
    for c in compounds {
        fp_values.extend_from_slice(&c.fingerprint);
    }
    let fp_values_array = Arc::new(UInt64Array::from(fp_values)) as Arc<dyn Array>;
    let fp_field = Arc::new(Field::new("item", DataType::UInt64, false));
    let fingerprint = FixedSizeListArray::try_new(fp_field, FP_WORDS as i32, fp_values_array, None)
        .map_err(|e| BitMakoError::Arrow(e))?;

    let batch = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(compound_ids),
            Arc::new(smiles),
            Arc::new(mw),
            Arc::new(logp),
            Arc::new(rot_bonds),
            Arc::new(heavy_atoms),
            Arc::new(ring_count),
            Arc::new(fingerprint),
        ],
    )
    .map_err(BitMakoError::Arrow)?;

    debug!("Built RecordBatch with {} rows", n);
    Ok(batch)
}

/// Write (or append to) a Lance dataset from an iterator of RecordBatches.
///
/// Uses lance's native Rust API. Must be called from an async context.
pub async fn write_lance_dataset(
    batches: Vec<RecordBatch>,
    output_path: &Path,
    append: bool,
) -> Result<()> {
    use lance::dataset::{Dataset, WriteMode, WriteParams};
    use arrow_array::RecordBatchIterator;

    if batches.is_empty() {
        return Ok(());
    }

    let schema = batches[0].schema();
    let mode = if append { WriteMode::Append } else { WriteMode::Create };

    let params = WriteParams {
        mode,
        max_rows_per_file: 1_000_000,
        max_rows_per_group: 1024,
        ..Default::default()
    };

    let batch_iter = RecordBatchIterator::new(
        batches.into_iter().map(Ok),
        schema,
    );

    Dataset::write(batch_iter, output_path.to_str().ok_or_else(|| {
        BitMakoError::Io(std::io::Error::new(std::io::ErrorKind::InvalidInput, "invalid path"))
    })?, Some(params))
    .await
    .map_err(|e| BitMakoError::Lance(e.to_string()))?;

    info!("Lance write complete: {}", output_path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etl::fingerprint::compute_morgan_fp;
    use crate::etl::properties::compute_properties;

    fn make_compound(smiles: &str, id: &str) -> ParsedCompound {
        ParsedCompound {
            compound_id: id.to_string(),
            smiles: smiles.to_string(),
            fingerprint: compute_morgan_fp(smiles),
            properties: compute_properties(smiles),
        }
    }

    #[test]
    fn test_build_record_batch() {
        let compounds = vec![
            make_compound("CCO", "Z001"),
            make_compound("c1ccccc1", "Z002"),
        ];
        let batch = build_record_batch(&compounds).unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.schema(), compound_schema());
    }

    #[test]
    fn test_schema_fields() {
        let schema = compound_schema();
        assert_eq!(schema.fields().len(), 8);
        assert_eq!(schema.field(0).name(), "compound_id");
        // fingerprint column should be FixedSizeList(UInt64, 16)
        match schema.field(7).data_type() {
            DataType::FixedSizeList(_, size) => assert_eq!(*size, 16),
            _ => panic!("Expected FixedSizeList"),
        }
    }
}
