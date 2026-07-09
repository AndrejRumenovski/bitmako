//! Full ETL pipeline integration test: builds a small bzip2-compressed TSV
//! fixture in-memory (no fixture file to commit — `bzip2` is already a main
//! dependency), runs it through the real `run_pipeline`, and checks the
//! resulting Lance dataset's row count, schema, and fingerprint agreement
//! with `compute_morgan_fp` called directly on the same SMILES.

use std::io::Write;

use bzip2::write::BzEncoder;
use bzip2::Compression;

use bitmako::etl::fingerprint::{compute_morgan_fp, FP_WORDS};
use bitmako::etl::writer::compound_schema;
use bitmako::etl::{run_pipeline, PipelineConfig};

/// (smiles, compound_id) rows for the fixture, matching Enamine REAL's
/// default column layout (smiles_col=0, id_col=1, tab-separated, header row).
const FIXTURE_ROWS: &[(&str, &str)] = &[
    ("CCO", "Z0000000001"),
    ("c1ccccc1", "Z0000000002"),
    ("CC(=O)Oc1ccccc1C(=O)O", "Z0000000003"),
    ("CNC(=O)c1ccccc1", "Z0000000004"),
    ("CCN(CC)CC", "Z0000000005"),
];

fn write_bz2_fixture() -> tempfile::NamedTempFile {
    let mut body = String::from("smiles\tid\n");
    for (smiles, id) in FIXTURE_ROWS {
        body.push_str(smiles);
        body.push('\t');
        body.push_str(id);
        body.push('\n');
    }

    let mut encoder = BzEncoder::new(Vec::new(), Compression::best());
    encoder.write_all(body.as_bytes()).unwrap();
    let compressed = encoder.finish().unwrap();

    let tmp = tempfile::Builder::new().suffix(".bz2").tempfile().unwrap();
    std::fs::write(tmp.path(), &compressed).unwrap();
    tmp
}

#[test]
fn run_pipeline_ingests_fixture_and_produces_a_matching_lance_dataset() {
    let fixture = write_bz2_fixture();
    let out_dir = tempfile::tempdir().unwrap();
    let dataset_path = out_dir.path().join("compounds.lance");

    let stats = run_pipeline(fixture.path(), &dataset_path, &PipelineConfig::default())
        .expect("pipeline must succeed on a well-formed fixture");

    assert_eq!(stats.total_lines, FIXTURE_ROWS.len());
    assert_eq!(stats.parsed_ok, FIXTURE_ROWS.len());
    assert_eq!(stats.parse_failures, 0);
    assert!(stats.batches_written >= 1);

    // Open the written Lance dataset and verify row count, schema, and that
    // every row's stored fingerprint matches compute_morgan_fp on the same
    // SMILES (proving the ETL's fingerprint column isn't silently diverging
    // from the fingerprinting code every other part of BitMako uses).
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    rt.block_on(async {
        use arrow_array::cast::AsArray;
        use arrow_array::types::UInt64Type;
        use arrow_array::{Array, FixedSizeListArray, StringArray};
        use futures::TryStreamExt;
        use lance::dataset::Dataset;

        let dataset = Dataset::open(dataset_path.to_str().unwrap()).await.unwrap();
        assert_eq!(dataset.schema().fields.len(), compound_schema().fields().len());

        let mut stream = dataset.scan().try_into_stream().await.unwrap();
        let mut seen_ids = Vec::new();
        while let Some(batch) = stream.try_next().await.unwrap() {
            assert_eq!(batch.num_columns(), 8);

            let ids = batch
                .column_by_name("compound_id")
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            let smiles_col = batch
                .column_by_name("smiles")
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            let fp_col = batch
                .column_by_name("fingerprint")
                .unwrap()
                .as_any()
                .downcast_ref::<FixedSizeListArray>()
                .unwrap();

            for row in 0..batch.num_rows() {
                seen_ids.push(ids.value(row).to_string());
                let smiles = smiles_col.value(row);
                let expected_fp = compute_morgan_fp(smiles);

                let values = fp_col.value(row);
                let u64_arr = values.as_primitive::<UInt64Type>();
                let mut stored_fp = [0u64; FP_WORDS];
                for (i, v) in u64_arr.values().iter().enumerate().take(FP_WORDS) {
                    stored_fp[i] = *v;
                }
                assert_eq!(
                    stored_fp, expected_fp,
                    "stored fingerprint for {smiles} diverges from compute_morgan_fp"
                );
            }
        }

        let mut expected_ids: Vec<String> = FIXTURE_ROWS.iter().map(|(_, id)| id.to_string()).collect();
        seen_ids.sort();
        expected_ids.sort();
        assert_eq!(seen_ids, expected_ids, "row count / identity mismatch after round-trip");
    });
}

#[test]
fn run_pipeline_counts_malformed_rows_as_parse_failures_not_a_hard_error() {
    // A blank-SMILES row is dropped by split_line/parsing rather than
    // aborting the whole ingest — the pipeline should keep going and report
    // it in parse_failures, not panic or return Err.
    let mut body = String::from("smiles\tid\n");
    body.push_str("CCO\tZ0000000001\n");
    body.push_str("\tZ0000000002\n"); // empty SMILES column
    body.push_str("c1ccccc1\tZ0000000003\n");

    let mut encoder = BzEncoder::new(Vec::new(), Compression::best());
    encoder.write_all(body.as_bytes()).unwrap();
    let compressed = encoder.finish().unwrap();
    let fixture = tempfile::Builder::new().suffix(".bz2").tempfile().unwrap();
    std::fs::write(fixture.path(), &compressed).unwrap();

    let out_dir = tempfile::tempdir().unwrap();
    let dataset_path = out_dir.path().join("compounds.lance");

    let stats = run_pipeline(fixture.path(), &dataset_path, &PipelineConfig::default()).unwrap();

    assert_eq!(stats.total_lines, 3);
    assert_eq!(stats.parsed_ok, 2);
    assert_eq!(stats.parse_failures, 1);
}
