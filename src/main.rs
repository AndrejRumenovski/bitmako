//! BitMako CLI entry point.
//!
//! Usage:
//!   bitmako ingest       --input <file.bz2> --output <dataset.lance>
//!   bitmako build-index  --lance <dataset.lance> --output <index.bitmako>
//!   bitmako search       --index <index.bitmako> --query <SMILES> --threshold 0.7 --top-k 100
//!   bitmako search       --index <index.bitmako> --query <SMILES> --mw-max 500 --logp-max 5

use std::path::PathBuf;
use std::process;

use tracing::{error, info};
use tracing_subscriber::{EnvFilter, fmt};

use bitmako::error::Result;
use bitmako::etl::{PipelineConfig, run_pipeline};
use bitmako::etl::fingerprint::{compute_morgan_fp, Fingerprint};
use bitmako::etl::reader::ReaderConfig;
use bitmako::index::builder::IndexBuilder;
use bitmako::index::IndexReader;
use bitmako::search::query::{PropertyField, PropertyFilter, SimilarityQuery};
use bitmako::search::Searcher;

fn init_tracing() {
    fmt()
        .with_env_filter(
            EnvFilter::from_default_env()
                .add_directive("bitmako=info".parse().unwrap()),
        )
        .with_target(false)
        .with_thread_ids(true)
        .json()
        .init();
}

fn main() {
    init_tracing();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: bitmako <ingest|build-index|search> [options]");
        process::exit(1);
    }

    let result = match args[1].as_str() {
        "ingest" => cmd_ingest(&args[2..]),
        "build-index" => cmd_build_index(&args[2..]),
        "search" => cmd_search(&args[2..]),
        other => {
            eprintln!("Unknown command: {}", other);
            process::exit(1);
        }
    };

    if let Err(e) = result {
        error!("Fatal error: {}", e);
        process::exit(1);
    }
}

/// `ingest --input <file.bz2> --output <dataset.lance> [--chunk-size N]`
fn cmd_ingest(args: &[String]) -> Result<()> {
    let input = PathBuf::from(require_flag(args, "--input"));
    let output = PathBuf::from(require_flag(args, "--output"));
    let chunk_size: usize = flag_value(args, "--chunk-size")
        .and_then(|s| s.parse().ok())
        .unwrap_or(100_000);

    let config = PipelineConfig {
        reader: ReaderConfig { chunk_size, ..ReaderConfig::default() },
        ..PipelineConfig::default()
    };

    info!("Starting ingestion: {:?} → {:?}", input, output);
    let stats = run_pipeline(&input, &output, &config)?;
    info!(
        "Ingestion complete: ok={} batches={} failures={}",
        stats.parsed_ok, stats.batches_written, stats.parse_failures
    );
    Ok(())
}

/// `build-index --lance <dataset.lance> --output <index.bitmako>`
fn cmd_build_index(args: &[String]) -> Result<()> {
    use arrow_array::{Array, cast::AsArray};
    use arrow_array::types::UInt64Type;

    let lance_path = PathBuf::from(require_flag(args, "--lance"));
    let index_path = PathBuf::from(require_flag(args, "--output"));

    info!("Building inverted index from Lance dataset: {:?}", lance_path);

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .map_err(bitmako::error::BitMakoError::Io)?;

    let mut builder = IndexBuilder::new();
    let mut doc_id: u32 = 0;

    rt.block_on(async {
        use lance::dataset::Dataset;
        use futures::TryStreamExt;

        let dataset = Dataset::open(
            lance_path.to_str().expect("non-UTF8 path"),
        )
        .await
        .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?;

        let mut stream = dataset
            .scan()
            .project(&["fingerprint"])
            .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?
            .try_into_stream()
            .await
            .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?;

        while let Some(batch) = stream
            .try_next()
            .await
            .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?
        {
            let fp_col = batch
                .column_by_name("fingerprint")
                .ok_or_else(|| bitmako::error::BitMakoError::IndexBuild(
                    "missing fingerprint column".into(),
                ))?;
            let list_arr = fp_col
                .as_any()
                .downcast_ref::<arrow_array::FixedSizeListArray>()
                .ok_or_else(|| bitmako::error::BitMakoError::IndexBuild(
                    "fingerprint not FixedSizeListArray".into(),
                ))?;

            for row in 0..list_arr.len() {
                let values = list_arr.value(row);
                let u64_arr = values.as_primitive::<UInt64Type>();
                let mut fp = [0u64; 16];
                for (i, v) in u64_arr.values().iter().enumerate().take(16) {
                    fp[i] = *v;
                }
                builder.add_compound(doc_id, &fp);
                doc_id += 1;
            }
        }
        Ok::<(), bitmako::error::BitMakoError>(())
    })?;

    let stats = builder.write_index(&index_path)?;
    info!(
        "Index built: compounds={} active_bits={} size_kb={}",
        stats.num_compounds,
        stats.non_empty_bits,
        stats.total_postings_bytes / 1024
    );
    Ok(())
}

/// `search --index <index.bitmako> --query <SMILES> --threshold 0.7 --top-k 50`
/// Optional: `--mw-max 500 --logp-max 5`
fn cmd_search(args: &[String]) -> Result<()> {
    let index_path = PathBuf::from(require_flag(args, "--index"));
    let query_smiles = require_flag(args, "--query");
    let threshold: f32 = flag_value(args, "--threshold")
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.7);
    let top_k: usize = flag_value(args, "--top-k")
        .and_then(|s| s.parse().ok())
        .unwrap_or(50);
    let mw_max: Option<f32> = flag_value(args, "--mw-max").and_then(|s| s.parse().ok());
    let logp_max: Option<f32> = flag_value(args, "--logp-max").and_then(|s| s.parse().ok());

    let index = IndexReader::open(&index_path)?;
    info!("Loaded index: {} compounds", index.num_compounds);
    info!("Query SMILES: '{}' threshold={} top_k={}", query_smiles, threshold, top_k);

    let query_fp = compute_morgan_fp(&query_smiles);
    let mut query = SimilarityQuery::new(query_fp, threshold, top_k);

    if let Some(max) = mw_max {
        query = query.with_filter(PropertyFilter {
            field: PropertyField::MolWeight,
            min: None,
            max: Some(max),
        });
    }
    if let Some(max) = logp_max {
        query = query.with_filter(PropertyFilter {
            field: PropertyField::LogP,
            min: None,
            max: Some(max),
        });
    }

    // NOTE: For a full integration, supply the fingerprint store from a
    // memory-mapped flat file written alongside the Lance dataset.
    let fp_store: Vec<Fingerprint> = Vec::new();
    let searcher = Searcher::open_from_index(index, fp_store);
    let results = searcher.search(&query)?;

    println!("Found {} results above threshold {}", results.len(), threshold);
    for (rank, (doc_id, score)) in results.iter().enumerate() {
        println!("  #{}: doc_id={} tanimoto={:.4}", rank + 1, doc_id, score);
    }
    Ok(())
}

fn require_flag(args: &[String], flag: &str) -> String {
    flag_value(args, flag).unwrap_or_else(|| {
        eprintln!("{} is required", flag);
        process::exit(1);
    })
}

fn flag_value(args: &[String], flag: &str) -> Option<String> {
    let pos = args.iter().position(|a| a == flag)?;
    args.get(pos + 1).cloned()
}
