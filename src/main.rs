//! BitMako CLI entry point.
//!
//! Usage:
//!   bitmako ingest         --input <file.bz2> --output <dataset.lance>
//!   bitmako build-index    --lance <dataset.lance> --output <index.bitmako>
//!   bitmako build-fp-store --lance <dataset.lance> --output <store.fp>
//!   bitmako search         --index <index.bitmako> --fp-store <store.fp> --query <SMILES> --threshold 0.7 --top-k 100
//!   bitmako search         --index <index.bitmako> --fp-store <store.fp> --query <SMILES> --mw-max 500 --logp-max 5

use std::path::PathBuf;
use std::process;

use tracing::{error, info};
use tracing_subscriber::{EnvFilter, fmt};

use bitmako::error::Result;
use bitmako::etl::{PipelineConfig, run_pipeline};
use bitmako::etl::fingerprint::compute_morgan_fp;
use bitmako::etl::reader::ReaderConfig;
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
        "build-fp-store" => cmd_build_fp_store(&args[2..]),
        "search" => cmd_search(&args[2..]),
        other => {
            eprintln!("Unknown command: {}", other);
            eprintln!("Usage: bitmako <ingest|build-index|build-fp-store|search> [options]");
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

/// `build-index --lance <dataset.lance> --output <index.bitmako> [--bits-per-pass N]`
///
/// Multi-pass streaming builder to avoid holding all 1024 posting lists in RAM.
/// Pass 0 collects compound popcounts (~1.4 GB for 1.4B compounds).
/// Subsequent passes each handle `bits_per_pass` bits, scanning Lance once per pass.
/// The offset table is written as a placeholder and patched after all passes complete.
///
/// Memory peak ≈ 1.4 GB (pops) + bits_per_pass × avg_list_bytes (~8 GB at default 32).
fn cmd_build_index(args: &[String]) -> Result<()> {
    use std::io::{Seek, SeekFrom, Write};
    use std::fs::OpenOptions;
    use arrow_array::{Array, cast::AsArray};
    use arrow_array::types::UInt64Type;
    use bitmako::etl::fingerprint::{fp_popcount, FP_BITS, FP_WORDS};
    use bitmako::index::posting_list::{PostingList, BLOCK_SIZE};

    let lance_path_str = require_flag(args, "--lance");
    let index_path = PathBuf::from(require_flag(args, "--output"));
    let bits_per_pass: usize = flag_value(args, "--bits-per-pass")
        .and_then(|s| s.parse().ok())
        .unwrap_or(32)
        .max(1)
        .min(FP_BITS);
    let num_passes = (FP_BITS + bits_per_pass - 1) / bits_per_pass;

    info!(
        "Building inverted index: {} bits/pass, {} total scans",
        bits_per_pass,
        num_passes + 1
    );

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .map_err(bitmako::error::BitMakoError::Io)?;

    // ---- Pass 0: collect one popcount byte per compound (~1.4 GB) ----
    info!("Pass 0: collecting compound popcounts...");
    let lance_str = lance_path_str.clone();
    let compound_pops: Vec<u8> = rt.block_on(async move {
        use lance::dataset::Dataset;
        use futures::TryStreamExt;

        let dataset = Dataset::open(&lance_str)
            .await
            .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?;
        let mut stream = dataset
            .scan()
            .project(&["fingerprint"])
            .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?
            .try_into_stream()
            .await
            .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?;

        let mut pops: Vec<u8> = Vec::new();
        while let Some(batch) = stream
            .try_next()
            .await
            .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?
        {
            let fp_col = batch
                .column_by_name("fingerprint")
                .ok_or_else(|| bitmako::error::BitMakoError::IndexBuild("missing fingerprint".into()))?;
            let list_arr = fp_col
                .as_any()
                .downcast_ref::<arrow_array::FixedSizeListArray>()
                .ok_or_else(|| bitmako::error::BitMakoError::IndexBuild("not FixedSizeListArray".into()))?;

            for row in 0..list_arr.len() {
                let values = list_arr.value(row);
                let u64_arr = values.as_primitive::<UInt64Type>();
                let mut fp = [0u64; FP_WORDS];
                for (i, v) in u64_arr.values().iter().enumerate().take(FP_WORDS) {
                    fp[i] = *v;
                }
                pops.push(fp_popcount(&fp) as u8);
            }
        }
        Ok::<Vec<u8>, bitmako::error::BitMakoError>(pops)
    })?;

    let num_compounds = compound_pops.len() as u32;
    info!("Pass 0 done: {} compounds", num_compounds);

    // ---- Write file header (placeholder offset table to be patched later) ----
    let mut out = OpenOptions::new()
        .read(true).write(true).create(true).truncate(true)
        .open(&index_path)?;

    out.write_all(b"BITMAKO1")?;
    out.write_all(&1u32.to_le_bytes())?;
    out.write_all(&num_compounds.to_le_bytes())?;
    out.write_all(&(FP_BITS as u32).to_le_bytes())?;
    // Byte 20: offset table (patched after all passes)
    const OFFSETS_POS: u64 = 20;
    out.write_all(&vec![0u8; FP_BITS * 8])?; // u64 per bit (v2 format)
    // Compound popcounts + 4-byte alignment padding
    out.write_all(&compound_pops)?;
    let pad = (4 - (num_compounds as usize % 4)) % 4;
    out.write_all(&vec![0u8; pad])?;

    // ---- Bit-group passes ----
    let mut bit_offsets = vec![0u64; FP_BITS];
    let mut posting_offset: u64 = 0;
    let mut total_bytes = 0usize;
    let mut non_empty_bits = 0usize;

    for pass in 0..num_passes {
        let bit_start = pass * bits_per_pass;
        let bit_end = (bit_start + bits_per_pass).min(FP_BITS);
        let n_bits = bit_end - bit_start;
        info!("Pass {}/{}: bits {}..{}", pass + 1, num_passes, bit_start, bit_end - 1);

        let lance_str = lance_path_str.clone();
        let doc_id_lists: Vec<Vec<u32>> = rt.block_on(async move {
            use lance::dataset::Dataset;
            use futures::TryStreamExt;

            let dataset = Dataset::open(&lance_str)
                .await
                .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?;
            let mut stream = dataset
                .scan()
                .project(&["fingerprint"])
                .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?
                .try_into_stream()
                .await
                .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?;

            let mut lists: Vec<Vec<u32>> = vec![Vec::new(); n_bits];
            let mut doc_id: u32 = 0;
            let word_start = bit_start / 64;
            let word_end = ((bit_end - 1) / 64) + 1;

            while let Some(batch) = stream
                .try_next()
                .await
                .map_err(|e| bitmako::error::BitMakoError::Lance(e.to_string()))?
            {
                let fp_col = batch
                    .column_by_name("fingerprint")
                    .ok_or_else(|| bitmako::error::BitMakoError::IndexBuild("missing fingerprint".into()))?;
                let list_arr = fp_col
                    .as_any()
                    .downcast_ref::<arrow_array::FixedSizeListArray>()
                    .ok_or_else(|| bitmako::error::BitMakoError::IndexBuild("not FixedSizeListArray".into()))?;

                for row in 0..list_arr.len() {
                    let values = list_arr.value(row);
                    let u64_arr = values.as_primitive::<UInt64Type>();
                    let words = u64_arr.values();

                    for word_idx in word_start..word_end.min(FP_WORDS) {
                        let mut bits = words.get(word_idx).copied().unwrap_or(0);
                        while bits != 0 {
                            let trailing = bits.trailing_zeros() as usize;
                            let global_bit = word_idx * 64 + trailing;
                            if global_bit >= bit_start && global_bit < bit_end {
                                lists[global_bit - bit_start].push(doc_id);
                            }
                            bits &= bits - 1;
                        }
                    }
                    doc_id += 1;
                }
            }
            Ok::<Vec<Vec<u32>>, bitmako::error::BitMakoError>(lists)
        })?;

        for (i, doc_ids) in doc_id_lists.into_iter().enumerate() {
            let bit = bit_start + i;
            bit_offsets[bit] = posting_offset;

            let num_blocks = (doc_ids.len() + BLOCK_SIZE - 1) / BLOCK_SIZE;
            let block_max_pop: Vec<u8> = (0..num_blocks)
                .map(|b| {
                    let s = b * BLOCK_SIZE;
                    let e = (s + BLOCK_SIZE).min(doc_ids.len());
                    doc_ids[s..e]
                        .iter()
                        .map(|&d| compound_pops.get(d as usize).copied().unwrap_or(0))
                        .max()
                        .unwrap_or(0)
                })
                .collect();

            if !doc_ids.is_empty() {
                non_empty_bits += 1;
            }
            let pl = PostingList { doc_ids, block_max_pop };
            let bytes = pl.serialize().map_err(bitmako::error::BitMakoError::Io)?;
            out.write_all(&bytes)?;
            posting_offset += bytes.len() as u64;
            total_bytes += bytes.len();
        }
    }

    // ---- Patch offset table now that all posting list positions are known ----
    out.seek(SeekFrom::Start(OFFSETS_POS))?;
    for &off in &bit_offsets {
        out.write_all(&off.to_le_bytes())?;
    }

    info!(
        "Index built: compounds={} active_bits={} size_kb={}",
        num_compounds, non_empty_bits, total_bytes / 1024
    );
    Ok(())
}

/// `build-fp-store --lance <dataset.lance> --output <store.fp>`
///
/// Writes every fingerprint to a flat binary file in Lance scan order, so the
/// resulting file can be memory-mapped and indexed directly by doc_id. Each
/// fingerprint is 16 little-endian u64 words (128 bytes). At 1.4B compounds the
/// output is ~174 GB — the scan order matches `build-index`, so doc_ids line up.
fn cmd_build_fp_store(args: &[String]) -> Result<()> {
    use std::io::{BufWriter, Write};
    use std::fs::File;
    use arrow_array::{Array, cast::AsArray};
    use arrow_array::types::UInt64Type;
    use bitmako::etl::fingerprint::FP_WORDS;

    let lance_path_str = require_flag(args, "--lance");
    let output = PathBuf::from(require_flag(args, "--output"));

    info!("Building flat fingerprint store: {} → {:?}", lance_path_str, output);

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .map_err(bitmako::error::BitMakoError::Io)?;

    let file = File::create(&output)?;
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, file);
    let mut count: u64 = 0;

    rt.block_on(async {
        use lance::dataset::Dataset;
        use futures::TryStreamExt;

        let dataset = Dataset::open(&lance_path_str)
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
                .ok_or_else(|| bitmako::error::BitMakoError::IndexBuild("missing fingerprint".into()))?;
            let list_arr = fp_col
                .as_any()
                .downcast_ref::<arrow_array::FixedSizeListArray>()
                .ok_or_else(|| bitmako::error::BitMakoError::IndexBuild("not FixedSizeListArray".into()))?;

            for row in 0..list_arr.len() {
                let values = list_arr.value(row);
                let u64_arr = values.as_primitive::<UInt64Type>();
                let words = u64_arr.values();

                let mut buf = [0u8; FP_WORDS * 8];
                for i in 0..FP_WORDS {
                    let w = words.get(i).copied().unwrap_or(0);
                    buf[i * 8..i * 8 + 8].copy_from_slice(&w.to_le_bytes());
                }
                writer.write_all(&buf).map_err(bitmako::error::BitMakoError::Io)?;
                count += 1;

                if count % 100_000_000 == 0 {
                    info!("Wrote {} million fingerprints", count / 1_000_000);
                }
            }
        }
        Ok::<(), bitmako::error::BitMakoError>(())
    })?;

    writer.flush().map_err(bitmako::error::BitMakoError::Io)?;
    info!(
        "Fingerprint store written: {} fingerprints ({} MB)",
        count,
        count * (FP_WORDS as u64 * 8) / (1024 * 1024)
    );
    Ok(())
}

/// `search --index <index.bitmako> --fp-store <store.fp> --query <SMILES> --threshold 0.7 --top-k 50`
/// Optional: `--mw-max 500 --logp-max 5`
fn cmd_search(args: &[String]) -> Result<()> {
    let index_path = PathBuf::from(require_flag(args, "--index"));
    let fp_store_path = PathBuf::from(require_flag(args, "--fp-store"));
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

    let fp_store = bitmako::search::fp_store::FpStore::open(&fp_store_path)?;
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
