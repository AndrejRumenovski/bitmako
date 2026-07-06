//! ETL pipeline: streaming bz2 → parallel parse → Arrow/Lance write.

pub mod fingerprint;
pub mod parser;
pub mod properties;
pub mod reader;
pub mod writer;

use std::path::Path;
use std::thread;

use crossbeam_channel::bounded;
use tracing::{info, span, Level};

use crate::error::{BitMakoError, Result};
use crate::etl::parser::parse_chunk_parallel;
use crate::etl::reader::{ReaderConfig, stream_bz2_file};
use crate::etl::writer::{build_record_batch, write_lance_dataset};

/// Full pipeline configuration
#[derive(Debug, Clone, Default)]
pub struct PipelineConfig {
    pub reader: ReaderConfig,
    /// Number of rayon threads (0 = use all CPUs)
    pub rayon_threads: usize,
}

/// Statistics collected during ingestion
#[derive(Debug, Default)]
pub struct IngestStats {
    pub total_lines: usize,
    pub parsed_ok: usize,
    pub parse_failures: usize,
    pub batches_written: usize,
}

/// Run the full ETL pipeline.
///
/// Streams the `.bz2` file, parses compound chunks in parallel with rayon,
/// builds Arrow RecordBatches, and writes Lance dataset fragments.
///
/// This function is synchronous at the top level (spawning an async runtime
/// internally for the Lance writes to avoid requiring callers to be async).
pub fn run_pipeline(
    input_path: &Path,
    output_path: &Path,
    config: &PipelineConfig,
) -> Result<IngestStats> {
    if config.rayon_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(config.rayon_threads)
            .build_global()
            .ok();
    }

    // Channel: streaming reader → parser workers
    // Backpressure capped to channel_capacity to bound memory usage
    let (tx, rx) = bounded::<Vec<crate::etl::reader::RawLine>>(config.reader.channel_capacity);

    let input_path_owned = input_path.to_path_buf();
    let reader_config = config.reader.clone();

    // Spawn streaming reader on a dedicated OS thread (blocking IO, not rayon)
    let reader_handle = thread::spawn(move || {
        stream_bz2_file(&input_path_owned, &reader_config, tx)
    });

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .map_err(BitMakoError::Io)?;

    let reader_config_2 = config.reader.clone();
    let mut stats = IngestStats::default();
    let mut append = false;

    // Receive chunks from reader, process in parallel, write to Lance
    for chunk in rx {
        let span = span!(Level::INFO, "parse_chunk", size = chunk.len());
        let _enter = span.enter();

        let (compounds, failed) = parse_chunk_parallel(chunk, &reader_config_2);
        stats.parsed_ok += compounds.len();
        stats.parse_failures += failed;

        if compounds.is_empty() {
            continue;
        }

        let batch = build_record_batch(&compounds)?;
        rt.block_on(write_lance_dataset(vec![batch], output_path, append))?;
        append = true;
        stats.batches_written += 1;

        info!(
            "Written batch {}: {} compounds ({} failed)",
            stats.batches_written, compounds.len(), failed
        );
    }

    stats.total_lines = reader_handle
        .join()
        .map_err(|_| BitMakoError::Io(std::io::Error::other("reader thread panicked")))??;

    info!(
        "Pipeline complete. total={} ok={} failed={} batches={}",
        stats.total_lines, stats.parsed_ok, stats.parse_failures, stats.batches_written
    );

    Ok(stats)
}
